# Dummy
**Maintainer**: *NAME* (*EMAIL*)
**Status**: Functional

### Description
This is a demonstrator module only, taking data every detector on the clipboard and plots the pixel hit positions.
It serves as template to create new modules.

### Parameters
No parameters are used from the configuration file.

### Plots produced
* Histogram of event numbers

For each detector the following plots are produced:

* 2D histogram of pixel hit positions

### Usage
```toml
[Dummy]

```
Parameters to be used in multiple modules can also be defined globally at the top of the configuration file. This is highly encouraged for parameters such as `DUT` and `reference`.