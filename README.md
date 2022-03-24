<p align="center">
  <a href="https://github.com/jlsteenwyk/orthosnap">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/orthosnap/master/docs/_static/img/logo.jpg" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/orthosnap/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/orthosnap/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/orthosnap/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://github.com/JLSteenwyk/orthosnap/actions" alt="Build">
            <img src="https://img.shields.io/github/workflow/status/jlsteenwyk/orthosnap/CI/master">
        </a>
        <a href="https://codecov.io/gh/JLSteenwyk/orthosnap">
          <img src="https://codecov.io/gh/JLSteenwyk/orthosnap/branch/master/graph/badge.svg?token=FX66FUET0L"/>
        </a>
        <a href="https://github.com/jlsteenwyk/orthosnap/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/orthosnap">
        </a>
        <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
            <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
                alt="follow on Twitter">
        </a>
        <br />
        <a href="https://pepy.tech/badge/orthosnap">
          <img src="https://static.pepy.tech/personalized-badge/orthosnap?period=total&units=international_system&left_color=grey&right_color=blue&left_text=PyPi%20Downloads">
        </a>
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/orthosnap/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/orthosnap">
        </a>
        <br />
        <a href="https://www.biorxiv.org/content/10.1101/2021.10.30.466607">
          <img src="https://zenodo.org/badge/DOI/10.1101/2021.10.30.466607.svg">
        </a>
        <a href="https://patreon.com/jlsteenwyk"><img src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Fshieldsio-patreon.vercel.app%2Fapi%3Fusername%3Djlsteenwyk%26type%3Dpatrons&style=flat" alt="Support me on Patreon" /></a>
    </p>
</p>

orthoSNAP: a tree splitting and pruning algorithm for retrieving single-copy orthologs from gene family trees.<br /><br />
If you found orthosnap useful, please cite *orthosnap: a tree splitting and pruning algorithm for retrieving single-copy orthologs from gene family trees*. Steenwyk et al. 2021, bioRxiv. doi: [10.1101/2021.10.30.466607](https://www.biorxiv.org/content/10.1101/2021.10.30.466607v1).
<br /><br />

---

<br />

This documentation covers downloading and installing orthosnap. Details about orthosnap usage including a tutorial are available on our [online documentation](https://jlsteenwyk.com/orthosnap/).

<br />

**Installation**

**If you are having trouble installing orthosnap, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [twitter](https://twitter.com/jlsteenwyk) to get help.**

To install using *pip*, we strongly recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:
```shell
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install orthosnap
pip install orthosnap
```
**Note, the virtual environment must be activated to use *orthosnap*.**

After using orthosnap, you may wish to deactivate your virtual environment and can do so using the following command:
```shell
# deactivate virtual environment
deactivate
```

<br />

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the following commands:
```shell
# download
git clone https://github.com/JLSteenwyk/orthosnap.git
cd orthosnap/
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install
make install
```
To deactivate your virtual environment, use the following command:
```shell
# deactivate virtual environment
deactivate
```
**Note, the virtual environment must be activated to use *orthosnap*.**

<br />

To install via anaconda, execute the follwoing command:

``` shell
conda install -c jlsteenwyk orthosnap
```
Visit here for more information: https://anaconda.org/jlsteenwyk/orthosnap


