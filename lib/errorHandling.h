#ifndef __errorHandling_h__
#define __errorHandling_h__

#include <iostream>
#include <stdexcept>
#include <string>

namespace error {

    /**
     * Processes exceptions and outputs appropriate error messages.
     * This function should be called inside a catch block.
     */
    inline void processException() {
        try {
            throw; // Re-throw the current exception to identify its type
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

}

#endif