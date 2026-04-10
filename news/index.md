# Changelog

## pc 0.3

## pc 0.2

#### new

- Document and maintain built-in case datasets for reproducible analysis
  ([\#61](https://github.com/stscl/pc/issues/61)).

- Extend `pc` generic with visualization capabilities
  ([\#55](https://github.com/stscl/pc/issues/55)).

- Provide `fnn` generic for *false nearest neighbors* method
  ([\#49](https://github.com/stscl/pc/issues/49)).

- Support linear trend removal for spatial cross-sectional data
  ([\#36](https://github.com/stscl/pc/issues/36)).

#### enhancements

- Support sampling with and without replacement in the bootstrapped
  version of pattern causality
  ([\#58](https://github.com/stscl/pc/issues/58)).

- Enable method dispatch compatibility in `pc` and `ops` generics via
  `...` ([\#43](https://github.com/stscl/pc/issues/43)).

- Improve handling of large-scale inputs
  ([\#33](https://github.com/stscl/pc/issues/33)).

#### bug fixes

- Fix unintended CCM-style xmap behavior in pattern causality estimation
  ([\#52](https://github.com/stscl/pc/issues/52)).

## pc 0.1

CRAN release: 2026-03-30

- First stable release ([\#24](https://github.com/stscl/pc/issues/24)).
