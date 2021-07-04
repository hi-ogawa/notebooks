#!/bin/bash

TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

jupyter nbconvert ${@} --to python --output-dir $TMP_DIR

touch $TMP_DIR/__init__.py
pylint $TMP_DIR
