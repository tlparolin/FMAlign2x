#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

class ThreadPool {
  public:
    explicit ThreadPool(size_t max_threads);
    ~ThreadPool();

    // Adds a new task to the pool
    void add_task(std::function<void()> task);

    // Waits for all tasks to finish
    void wait_for_tasks();

    // Waits for all tasks to finish and finalizes the pool
    void shutdown();

  private:
    void worker_thread();

    size_t max_threads_;
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable cv_done_;
    bool quitting_ = false;
    size_t active_tasks_ = 0;
};

#endif // THREAD_POOL_H
