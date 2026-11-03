# Builds ue (sketches, statistics, theory), calibrate (C and b of the schemes, from simulated
# copies) and copies (suffix array and LCP of a genome).
# copies needs libdivsufsort; the first `make` downloads and installs it into
# src/third_party/libdivsufsort (static, 32- and 64-bit index), which needs git or curl,
# cmake and a C compiler.
#
#   make -C src              ue, calibrate and copies
#   make -C src ue calibrate ue and calibrate only (no libdivsufsort needed)
#   make -C src divsufsort   libdivsufsort only
#   make -C src clean        remove ue, calibrate and copies (keeps libdivsufsort)
#   make -C src distclean    also remove third_party/

CXXFLAGS ?= -O3 -std=c++17 -Wall
FPFLAGS  := -ffp-contract=off   # no fused multiply-add: ue's floating point matches numpy
JOBS     ?= 4

DS_PREFIX  := $(CURDIR)/third_party/libdivsufsort
DS_URL     := https://github.com/y-256/libdivsufsort.git
DS_TARBALL := https://github.com/y-256/libdivsufsort/archive/refs/heads/master.tar.gz
DS_LIB     := $(DS_PREFIX)/lib/libdivsufsort64.a

.PHONY: all divsufsort clean distclean
all: ue calibrate copies

ue: main.cpp sketches.cpp theory.cpp sketches.h theory.h
	$(CXX) $(CXXFLAGS) $(FPFLAGS) -o $@ main.cpp sketches.cpp theory.cpp -lz

calibrate: calibrate.cpp sketches.cpp theory.cpp sketches.h theory.h
	$(CXX) $(CXXFLAGS) $(FPFLAGS) -o $@ calibrate.cpp sketches.cpp theory.cpp

copies: copies.cpp $(DS_LIB)
	$(CXX) $(CXXFLAGS) -I$(DS_PREFIX)/include -o $@ copies.cpp -L$(DS_PREFIX)/lib -ldivsufsort64 -lz

divsufsort: $(DS_LIB)

# CMAKE_POLICY_VERSION_MINIMUM lets CMake 4 configure libdivsufsort's old CMakeLists (older CMake ignores it)
$(DS_LIB):
	@set -e; \
	command -v cmake >/dev/null || { echo "error: libdivsufsort needs cmake (e.g. apt install cmake, brew install cmake, module load cmake)"; exit 1; }; \
	w=$$(mktemp -d); \
	echo ">> fetching libdivsufsort"; \
	git clone --depth 1 $(DS_URL) $$w/src >/dev/null 2>&1 || { \
		mkdir -p $$w/src; curl -fsSL $(DS_TARBALL) | tar xz -C $$w/src --strip-components=1; }; \
	echo ">> building libdivsufsort (static, 32- and 64-bit) into $(DS_PREFIX)"; \
	mkdir -p $$w/build; cd $$w/build; \
	cmake -Wno-dev -Wno-deprecated --no-warn-unused-cli -DCMAKE_BUILD_TYPE=Release -DBUILD_DIVSUFSORT64=ON -DBUILD_SHARED_LIBS=OFF -DBUILD_EXAMPLES=OFF \
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_INSTALL_PREFIX=$(DS_PREFIX) $$w/src >/dev/null; \
	$(MAKE) -j $(JOBS) >/dev/null; \
	$(MAKE) install >/dev/null; \
	cd $(CURDIR); rm -rf $$w; \
	test -f $@ && echo ">> libdivsufsort installed"

clean:
	rm -f ue calibrate copies

distclean: clean
	rm -rf third_party
