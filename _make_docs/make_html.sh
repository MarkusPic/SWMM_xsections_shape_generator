#!/usr/bin/env bash

make html
mv ../docs/html ../docs_
rm -r ../docs
mv ../docs_ ../docs