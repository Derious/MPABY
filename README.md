# MPABY: Faster Secure Multi-Party Computation using Mixed-Mode Framework

## Specifications

- OS: Linux x64
- Language: C++
- Requires: emp-tools,  emp-ot, Eigen-3.4.0

## Compiling

To compile and test the conversion protocol, into the root dir and do the following: 

```
  $ mkdir build && cd build
  $ cmake ../
  $ make test_test_A2B_sh
  $ make test_test_Bit2A_sh
```

## Running

If you want to test the code in local machine, into the build dir and  type

```
  $../run ./bin/test_test_A2B_sh //when nP=2 which means it involves 2 parties
```

If you want to test the code over two machine, into the build dir and type

```
  $./bin/[binaries] 1 12345  //on one machine and

  $./bin/[binaries] 2 12345  //on the other.
```

You can modify the IP in the **MPABY/util_cmpc_config.h** to communicate with multiple machines
