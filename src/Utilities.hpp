#pragma once
#include <cstdio>
#include <string>
#include <chrono>
#include "Params.hpp"
#include <string>
#include <cstdarg>

namespace Parallel {
    inline int nRanks   = -1;
    inline int iRank    = -1;
    inline int nThreads = -1;

    inline bool isRoot() { return iRank == 0; }
}

//============================================================================

class Logger {
private:
    void applyIndentDelta(int delta) {
        indentLevel += delta;
        if (indentLevel < 0) indentLevel = 0;
    }

    std::string getIndent() const {
        std::string indent;
        indent.reserve(indentLevel * indentString.size());
        for (int i = 0; i < indentLevel; ++i) {
            indent += indentString;
        }
        return indent;
    }

    template <typename... Args>
    std::string formatToString(const char* fmt, Args... args) const {
        if constexpr (sizeof...(Args) == 0) {
            return std::string(fmt);
        } else {
            int size = std::snprintf(nullptr, 0, fmt, args...);
            if (size <= 0) return std::string(fmt);

            std::string buf(size, '\0');
            std::snprintf(buf.data(), size + 1, fmt, args...);
            return buf;
        }
    }

public:
    int indentLevel = 0;
    std::string indentString = "     ";

    template <typename... Args>
    Logger& print(int indentDelta, const char* fmt, Args... args) {
        if (!Parallel::isRoot()) return *this;

        if (indentDelta < 0) applyIndentDelta(indentDelta);

        if constexpr (sizeof...(Args) == 0) {
            std::printf("%s", fmt);
        } else {
            std::printf(fmt, args...);
        }

        std::fflush(stdout);

        if (indentDelta > 0) applyIndentDelta(indentDelta);

        return *this;
    }

    template <typename... Args>
    Logger& print(const char* fmt, Args... args) {
        return print(0, fmt, args...);
    }

    template <typename... Args>
    Logger& formPrint(int indentDelta, const char* fmt, Args... args) {
        if (!Parallel::isRoot()) return *this;

        if (indentDelta < 0) applyIndentDelta(indentDelta);

        std::string msg = formatToString(fmt, args...);

        bool hasNewline = !msg.empty() && msg.back() == '\n';
        if (hasNewline) {
            msg.pop_back();
        }

        std::string indent = getIndent();
        int indentLen = static_cast<int>(indent.size());
        int msgLen = static_cast<int>(msg.size());
        int pad = Params::w - (indentLen + msgLen);

        std::printf("%s%s", indent.c_str(), msg.c_str());
        if (pad > 0) {
            std::printf("%*s", pad, "");
        }

        // 5. Re-attach newline if present
        if (hasNewline) {
            std::printf("\n");
        }

        std::fflush(stdout);

        if (indentDelta > 0) applyIndentDelta(indentDelta);

        return *this;
    }

    template <typename... Args>
    Logger& formPrint(const char* fmt, Args... args) {
        return formPrint(0, fmt, args...);
    }
};

inline Logger logg;

//============================================================================

class Timer {
private:
    static inline std::chrono::high_resolution_clock::time_point startTime;

public:
    static void start(const char* msg = nullptr) {
        if (msg) {
            logg.formPrint(msg);
        }
        startTime = std::chrono::high_resolution_clock::now();
    }

    static void end() {
        if (!Parallel::isRoot()) return;

        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = endTime - startTime;
        
        logg.print("done. -- %8.4f s\n", elapsed.count());
    }

};