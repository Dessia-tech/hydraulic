#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:35:47 2021

@author: steven
"""

import json

MIN_FILE_COVERAGE = 79.0
MIN_MODULE_COVERAGE = 71.0

untracked_modules = [

]

print("untracked modules:", untracked_modules)

with open("coverage.json", "r") as file:
    d = json.load(file)

print("total covered", d["totals"]["percent_covered"], "%")
assert d["totals"]["percent_covered"] >= MIN_MODULE_COVERAGE

for file_name, data in d["files"].items():
    print(file_name, data["summary"]["percent_covered"], "%")
    if "/".join(file_name.split("/")[-2:]) in untracked_modules:
        print(file_name, "-> in untrack list")
    else:
        # print('Testing if {} is above {}'.format(file_name, MIN_FILE_COVERAGE))
        assert (data["summary"]["percent_covered"]) >= MIN_FILE_COVERAGE
