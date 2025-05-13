#pragma once
#if defined(_WIN32)
    #ifdef ctrecon_math_EXPORTS
        #define ctrecon_math_API __declspec(dllexport)
    #else
        #define ctrecon_math_API __declspec(dllimport)
    #endif
#else
    #ifdef ctrecon_math_EXPORTS
        #define ctrecon_math_API __attribute__((visibility("default")))
    #else
        #define ctrecon_math_API
    #endif
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    ctrecon_math_API const char* get_version();
    ctrecon_math_API int test();

#ifdef __cplusplus
}
#endif

#undef ctrecon_math_API
