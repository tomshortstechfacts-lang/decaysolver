// T1 — l'en-tête de provenance doit être présent, commenté (`#`) et horodaté en UTC ISO 8601.

#include <decaysolver/provenance.hpp>
#include <decaysolver/version.hpp>

#include <catch2/catch_test_macros.hpp>

#include <regex>
#include <sstream>
#include <string>

TEST_CASE("provenance: horodatage UTC ISO 8601", "[provenance][T1]") {
    const std::regex iso8601_utc(R"(^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$)");
    REQUIRE(std::regex_match(decaysolver::utc_timestamp_iso8601(), iso8601_utc));
}

TEST_CASE("provenance: en-tete des sorties", "[provenance][T1]") {
    const std::string header = decaysolver::provenance_header();

    SECTION("chaque ligne est un commentaire CSV") {
        std::istringstream lines(header);
        std::string line;
        int count = 0;
        while (std::getline(lines, line)) {
            REQUIRE(!line.empty());
            REQUIRE(line.front() == '#');
            ++count;
        }
        REQUIRE(count == 5);
    }
    SECTION("contient la version, le SHA git et le compilateur") {
        REQUIRE(header.find(decaysolver::version::string) != std::string::npos);
        REQUIRE(header.find("git_sha: ") != std::string::npos);
        REQUIRE(header.find(decaysolver::version::compiler) != std::string::npos);
    }
}
