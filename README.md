# CMIMSUNIT

## MIMSUNIT API

See `cmims/src/mims_unit.h` for function API, and `cmims/tst/demo.c` for usage examples.

---

## How to run demo and tests

### Create folder for test files

```bash
mkdir data && cd data
```

Download `mims_unit.zip` (ask Aryton) and unzip in `cmims/data`

### Compile and run tests

```bash
cd cmims/tst && make clean
./demo
```

#### With Cmake

Run from the root directory, the build system will be in `build` folder.

```bash
# Generate the build system
cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -S . -B ./build -T host=x64 -A x64
```

Then compile using the build system, the compiled binary will be in `build/Debug` folder.

```bash
cmake --build ./build --config Debug --target ALL_BUILD -j 14 --
```

To compile for release, replace the `Debug` word in the above two commands with `Release`, the compiled binary will be in `build/Release` folder.

### Run mims_unit with different inputs

Modify `main()` in `tst/demo.c`
