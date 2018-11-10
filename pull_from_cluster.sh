#!/bin/bash
rsync -av --exclude-from 'rsync-exclude.txt' cluster1:quasi-laplace-real-analysis/ ./
