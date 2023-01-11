## MIMSUNIT API

See `cmims/src/mims_unit.h` for function API, and `cmims/tst/demo.c` for usage examples.

---

## How to run demo and tests

### Create folder for test files

```
mkdir data && cd data
```

Download `mims_unit.zip` (ask Aryton) and unzip in `cmims/data`

### Compile and run tests

```
cd cmims/tst && make clean
./demo
```

### Run mims_unit with different inputs

Modify `main()` in `tst/demo.c`
