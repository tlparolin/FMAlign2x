#include "thread_pool.h"

ThreadPool::ThreadPool(size_t max_threads) : max_threads_(max_threads) {
    for (size_t i = 0; i < max_threads_; ++i) {
        workers_.emplace_back(&ThreadPool::worker_thread, this);
    }
}

ThreadPool::~ThreadPool() { shutdown(); }

void ThreadPool::add_task(std::function<void()> task) {
    {
        std::lock_guard lock(mutex_);
        tasks_.push(std::move(task));
    }
    cv_.notify_one();
}

void ThreadPool::shutdown() {
    quitting_ = true;
    cv_.notify_all();

    for (auto &worker : workers_) {
        if (worker.joinable())
            worker.join();
    }
    workers_.clear();
}

void ThreadPool::worker_thread() {
    while (true) {
        std::function<void()> task;

        {
            std::unique_lock lock(mutex_);
            cv_.wait(lock, [this] { return quitting_ || !tasks_.empty(); });

            if (quitting_ && tasks_.empty())
                return;

            task = std::move(tasks_.front());
            tasks_.pop();
        }

        task();
    }
}
