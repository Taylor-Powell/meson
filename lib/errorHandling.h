#ifndef __errorHandling_h__
#define __errorHandling_h__

#include <iostream>
#include <stdexcept>
#include <future>
#include <system_error>
#include <exception>
#include <string>

namespace error {
    /**
     * Processes exceptions and outputs appropriate error messages.
     * This function should be called inside a catch block.
     * Pulled from The C++ Standard Library 2nd edition (page 50-51)
     */
    inline void processException() {
        try {
            throw; // Re-throw the current exception to identify its type
        }
        catch (const std::ios_base::failure& e) {
            std::cerr << "\n\nI/O EXCEPTION:\n" << e.what() << std::endl;
            processCodeException(e);
        }
        catch (const std::system_error& e) {
            std::cerr << "\n\nSystem EXCEPTION:\n" << e.what() << std::endl;
            processCodeException(e);
        }
        catch (const std::future_error& e) {
            std::cerr << "\n\nFuture EXCEPTION:\n" << e.what() << std::endl;
            processCodeException(e);
        }
        catch (const std::bad_alloc& e) {
            std::cerr << "\n\nBAD ALLOC EXCEPTION:\n" << e.what() << std::endl;
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "\n\nInvalid Argument Error:\n" << e.what() << std::endl;
        }
        catch (const std::domain_error& e) {
            std::cerr << "\n\nDomain Error:\n" << e.what() << std::endl;
        }
        catch (const std::runtime_error& e) {
            std::cerr << "\n\nRuntime Error:\n" << e.what() << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "\n\nGeneral Exception:\n" << e.what() << std::endl;
        }
        catch (const std::string& e) {
            std::cerr << "\n\nString Error:\n" << e << std::endl;
        }
        catch (...) {
            std::cerr << "\n\nUnknown Error occurred.\n";
        }
    }

    template <typename T>
    void processCodeException(const T& e) {
        auto c = e.code();
        std::cerr << "- category: " << c.category().name() << std::endl;
        std::cerr << "-value: " << c.value() << std::endl;
        std::cerr << "- message: " << c.message() << std::endl;
        std::cerr << "- def category: " << c.default_error_condition().name() << std::endl;
        std::cerr << "- def value: " << c.default_error_condition().value() << std::endl;
        std::cerr << "- def message: " << c.default_error_condition().message() << std::endl;
    }
}

#endif