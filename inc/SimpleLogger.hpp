/**
 * @file SimpleLogger.hpp
 * @author your name (you@domain.com)
 * @brief a simple logger that can be included and used directly
 * @version 0.1
 * @date 2022-08-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef SIMPLELOGGER_HPP
#define SIMPLELOGGER_HPP

#include <iostream>
#include <string>
#include <vector>

class SimpleLogger {
    public:
    int errSeverity = 9;

    struct logMsg{
        int severity;
        std::string msg;
    };

    // use 0 as nothing to 9 as everything / severity level (high big number, low small number)
    // 9 is fails, 0 ist just some status and 4 or so is x result
    int logLevel = 0;
    bool printLogs = true;
    std::vector<logMsg> logList;

    void log(int severity, std::string msg) {
        struct logMsg lm = {severity, msg};
        this->logList.push_back(lm);
        if (printLogs && severity >= this->logLevel) {
            std::cout << "[" << severity << "]: " << msg << std::endl;
        }
    }

    void addErrLog(std::string text) {
        log(errSeverity, text);
    }
    std::vector<logMsg> getErrLogs(int severity=-1) {
        severity = severity > 0 ? severity: errSeverity;
        std::vector<logMsg> errLogs;
        for (auto msg: this->logList) {
            if (msg.severity >= severity) {
                errLogs.push_back(msg);
            }
        }
        return errLogs;
    }
    void printErrLogs(int n=-1) {
        auto errLogs = getErrLogs();
        if (n == -1) {
            for (auto msg: errLogs) {
                std::cout << msg.msg << std::endl;
                return;
            }
        }
        if (n < 0) {
            // for negative numbers != -1 go from the end
            n += errLogs.size();
            // if n is still negative make it invalid
            n = (n < 0) ? errLogs.size()+1 : errLogs.size();
        }
        // at this point n is > 0
        if (static_cast<uint>(n) >= errLogs.size()) {
            cout << "invalid error idx. idx: " << n << " number of errors is: " << errLogs.size() << endl;
            return;
        }
        std::cout << "error (" << n << "): " << errLogs.at(n).msg << std::endl;
    }
    void clearErrLogs(int severity=-1) {
        severity = severity > 0 ? severity: errSeverity;
        this->logList.erase(std::remove_if(this->logList.begin(), this->logList.end(), [severity](logMsg msg){return msg.severity > severity;}), this->logList.end());
        // https://www.cplusplus.com/reference/vector/vector/clear/
        // reallocation is not guaranteed to happen, and the vector capacity is not guaranteed to change due to ca
        this->logList.shrink_to_fit();
    }
    void printLastErr() {
        printErrLogs(this->getErrLogs().size() - 1);
    }
};
// obv its always inline but this way its easyer readayble maybe
// usage is include ...
// logger.log()
// inline ensures that its always the same object / pointer
inline SimpleLogger logger;

#endif // SIMPLELOGGER_HPP