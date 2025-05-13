#pragma once
#include <array>
#include <atomic>
#include <future>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <vector>

struct irender
{
    virtual void render() {}
};

class debugger
{
    std::atomic<std::shared_ptr<irender>> render;

public:
    debugger(std::shared_ptr<irender> render) : render(render)
    {
        if (render == nullptr)
            throw std::runtime_error("render is nullptr");
    }
    ~debugger() = default;
    void initialize();
    void destory();
    int execute();

public:
    void* window = nullptr;
    std::future<void> window_main;
    std::atomic_flag window_inited = ATOMIC_FLAG_INIT;
    std::atomic_flag window_running = ATOMIC_FLAG_INIT;
    std::atomic_flag window_destoryed = ATOMIC_FLAG_INIT;
};
