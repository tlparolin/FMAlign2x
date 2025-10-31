#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

class ThreadPool {
  public:
    explicit ThreadPool(size_t max_threads);
    ~ThreadPool();

    void add_task(std::function<void()> task);
    void wait_for_tasks();
    void shutdown();

  private:
    void worker_thread(size_t index);
    std::optional<std::function<void()>> steal_task(size_t thief_index);

    size_t max_threads_;
    std::vector<std::thread> workers_;

    // Cada thread tem seu deque de tarefas
    std::vector<std::deque<std::function<void()>>> task_queues_;
    // Mutex individual para cada fila (protege acesso ao deque desta thread)
    std::vector<std::mutex> queue_mutexes_;

    std::condition_variable_any cv_; // Usar cv_any pois teremos vários mutexes
    std::atomic<bool> quitting_{false};

    std::atomic<size_t> active_tasks_{0};
    std::mutex wait_mutex_;
    std::condition_variable cv_done_;
};

#endif // THREAD_POOL_H