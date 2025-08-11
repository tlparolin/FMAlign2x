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

    // Adiciona uma tarefa para execução
    void add_task(std::function<void()> task);

    // Espera todas as tarefas finalizarem e destrói o pool
    void shutdown();

  private:
    void worker_thread();

    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;

    std::mutex mutex_;
    std::condition_variable cv_;
    std::atomic<bool> quitting_{false};
    size_t max_threads_;
};

#endif // THREAD_POOL_H
