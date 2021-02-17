# alertstack

[![Build Status](https://travis-ci.org/icecube/alertstack.svg?branch=master)](https://travis-ci.org/icecube/alertstack) [![Coverage Status](https://coveralls.io/repos/github/icecube/alertstack/badge.svg)](https://coveralls.io/github/icecube/alertstack)  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/icecube/alertstack/HEAD)

Code to test correlations of catalogues and neutrino alerts, created by [@robertdstein](https://github.com/robertdstein) and [@clagunas](https://github.com/clagunas).

*This code is still under development, so is not endorsed by the IceCube collaboration for external analysis.*

## Get Started

Try running *alertstack* at:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/icecube/alertstack/HEAD)

and then look at examples to get an idea.

## Install locally

To install locally:

`git clone https://github.com/icecube/alertstack.git`
`pip install -e alertstack/`

Pip installation with pypi will soon be supported too.

## Use full skymaps

*alertstack* runs both with published neutrino information or internal likelihood skymaps. Icecube Collaboration members can run:

`export NU_SKYMAP_DIR=/path/to/healpix/files `

to import and use Healpix files. 