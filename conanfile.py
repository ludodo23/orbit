import os
import re
from pathlib import Path
from conan import ConanFile
from conan.tools.files import copy
from conan.tools.cmake import cmake_layout, CMake, CMakeToolchain, CMakeDeps
from conan.tools.build import check_min_cppstd


class OrbitConan(ConanFile):
    name = "orbit"
    settings = "os", "arch", "compiler", "build_type"
    exports_sources = "include/*", "tests/*"
    no_copy_source = True
    package_type = "header-library"

    def set_version(self):
        header = Path(self.recipe_folder) / "include" / "orbital_elements.hpp"
        content = header.read_text()
        version = re.search(r'ORBIT_VERSION\s+"(.+)"', content)
        if version:
            self.version = version.group(1)

    def build_requirements(self):
        self.test_requires("catch2/3.4.0")

    def validate(self):
        check_min_cppstd(self, 14)

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        if self.conf.get("tools.build:skip_test", default=False):
            return
        cmake = CMake(self)
        cmake.configure()   # prend la racine, plus build_script_folder
        cmake.build()
        cmake.ctest(cli_args=["--output-on-failure"])

    def package(self):
        copy(self, "*.hpp", self.source_folder, self.package_folder)

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []

    def package_id(self):
        self.info.clear()