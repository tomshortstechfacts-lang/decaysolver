#include <decaysolver/provenance.hpp>
#include <decaysolver/version.hpp>

#include <array>
#include <ctime>
#include <sstream>
#include <string>

namespace decaysolver {

std::string utc_timestamp_iso8601() {
    const std::time_t now = std::time(nullptr);
    std::tm utc{};
#if defined(_WIN32)
    gmtime_s(&utc, &now);
#else
    gmtime_r(&now, &utc);
#endif
    std::array<char, 32> buffer{};
    const std::size_t written =
        std::strftime(buffer.data(), buffer.size(), "%Y-%m-%dT%H:%M:%SZ", &utc);
    return std::string(buffer.data(), written);
}

std::string provenance_header() {
    std::ostringstream out;
    out << "# decaysolver: " << version::string << '\n';
    out << "# git_sha: " << version::git_sha << (version::git_dirty ? " (arbre modifié)" : "")
        << " (capturé à la configuration CMake)\n";
    out << "# compiler: " << version::compiler << '\n';
    out << "# build_type: " << version::build_type << '\n';
    out << "# generated_utc: " << utc_timestamp_iso8601() << '\n';
    return out.str();
}

} // namespace decaysolver
