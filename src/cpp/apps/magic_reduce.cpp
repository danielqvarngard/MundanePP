#include <iostream>
#include <string_view>

namespace {

void print_help() {
    std::cout
        << "magic_reduce - Clifford and magic circuit preconditioner\n\n"
        << "Usage:\n"
        << "  magic_reduce --help\n\n"
        << "The optimization algorithm is not implemented yet.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc == 1 || (argc == 2 && std::string_view{argv[1]} == "--help")) {
        print_help();
        return 0;
    }

    std::cerr << "error: unsupported arguments; use --help\n";
    return 1;
}
