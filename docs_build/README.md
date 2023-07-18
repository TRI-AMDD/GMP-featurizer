For documentation generation, start with `sphinx-quickstart`, then modify make.bat and the Makefile to target `../docs`.

For apidoc generation, ensure the `sphinx.ext.napoleon` extension is enabled in `source/conf.py`. Then run ` sphinx-apidoc --separate -d 7 -o source/ ../GMPFeaturizer`.

Finally, for generating the HTML, run `make html` in this folder, which should produce the docs, then push to do the docs branch, which triggers the github build.