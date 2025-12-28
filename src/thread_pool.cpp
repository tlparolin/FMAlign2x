#include "thread_pool.h"

ThreadPool::ThreadPool(size_t max_threads) : max_threads_(max_threads), task_queues_(max_threads), queue_mutexes_(max_threads) {
    for (size_t i = 0; i < max_threads_; ++i) {
        workers_.emplace_back(&ThreadPool::worker_thread, this, i);
    }
}

ThreadPool::~ThreadPool() { shutdown(); }

void ThreadPool::add_task(std::function<void()> task) {
    // Distributes tasks to thread queues in a simple round-robin fashion.
    static std::atomic<size_t> next_queue{0};
    size_t queue_index = next_queue++ % max_threads_;

    {
        std::lock_guard lock(queue_mutexes_[queue_index]);
        task_queues_[queue_index].push_front(std::move(task)); // Push front to improve stealing
    }
    cv_.notify_all();
}

void ThreadPool::wait_for_tasks() {
    std::unique_lock lock(wait_mutex_);
    cv_done_.wait(lock, [this] {
        // Wait until no tasks and no tasks in progress.
        bool empty = true;
        for (size_t i = 0; i < max_threads_; ++i) {
            std::lock_guard qlock(queue_mutexes_[i]);
            if (!task_queues_[i].empty()) {
                empty = false;
                break;
            }
        }
        return empty && active_tasks_ == 0;
    });
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

std::optional<std::function<void()>> ThreadPool::steal_task(size_t thief_index) {
    // Try stealing from the queue of other threads starting with next.
    for (size_t i = 0; i < max_threads_; ++i) {
        size_t victim = (thief_index + i + 1) % max_threads_; // avoids stealing from the queue itself

        std::lock_guard lock(queue_mutexes_[victim]);
        if (!task_queues_[victim].empty()) {
            auto task = std::move(task_queues_[victim].back());
            task_queues_[victim].pop_back();
            return task;
        }
    }
    return {}; // nothing to stealing
}

void ThreadPool::worker_thread(size_t index) {
    while (true) {
        std::function<void()> task;

        // First, try searching within your own queue.
        {
            std::unique_lock lock(queue_mutexes_[index]);
            if (task_queues_[index].empty()) {
                // Try to steal from another queue
                lock.unlock();

                auto stolen = steal_task(index);
                if (stolen.has_value()) {
                    task = std::move(stolen.value());
                }
            } else {
                task = std::move(task_queues_[index].front());
                task_queues_[index].pop_front();
            }
        }

        if (!task && quitting_) {
            return; // Leave if it's almost finished and there are no tasks left.
        }

        if (!task) {
            // Nothing to do now, wait for notification.
            std::unique_lock<std::mutex> wait_lock(wait_mutex_);
            cv_.wait_for(wait_lock, std::chrono::milliseconds(10));
            continue;
        }

        active_tasks_++;
        task();
        active_tasks_--;

        if (active_tasks_ == 0) {
            std::lock_guard<std::mutex> wait_lock(wait_mutex_);
            cv_done_.notify_all();
        }
    }
}
