#!/usr/bin/python

"""
Originally from here: https://gist.github.com/daler/5473540
- Submit example data to REVIGO server (http://revigo.irb.hr/)
- Download and run R script for creating the treemap
"""

import os
import sys
import urllib
import mechanize

url = "http://revigo.irb.hr/"

# RobustFactory because REVIGO forms not well-formatted
br = mechanize.Browser()

# Read in data
# txt = open('Fr4VsInput_posLFC_revigo_in.txt').read()
txt = open(sys.argv[1]).read()
label = sys.argv[2]

# Encode and request
data = {'inputGoList': txt}
br.open(url, data=urllib.urlencode(data))

# Submit form
br.select_form(name="submitToRevigo")
response = br.submit()


for link in br.links(url_regex='toR_treemap.jsp'):
    if "toR_treemap.jsp?table=1" in link.url:
        br.follow_link(url="toR_treemap.jsp?table=1")
        rscript = label + '_treemap_bp.R'
        with open(rscript, 'w') as f:
            f.write(br.response().read())
    if "toR_treemap.jsp?table=2" in link.url:
        if os.path.exists(label + '_treemap_bp.R'):
            br.back()
        br.follow_link(url="toR_treemap.jsp?table=2")
        rscript = label + '_treemap_cc.R'
        with open(rscript, 'w') as f:
            f.write(br.response().read())
    if "toR_treemap.jsp?table=3" in link.url:
        if os.path.exists(label + '_treemap_cc.R'):
            br.back()
        br.follow_link(url="toR_treemap.jsp?table=3")
        rscript = label + '_treemap_mf.R'
        with open(rscript, 'w') as f:
            f.write(br.response().read())
