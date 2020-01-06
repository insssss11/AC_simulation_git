#!/usr/bin/env python
# -*- coding: utf-8 -*-

from subprocess import call
desDir = './result'
srcDir = 'simpson9794@nuclear.korea.ac.kr:/home/simpson9794/AC_simulation_git/result/'
call(['rsync', '-ravzP', '--delete', srcDir, desDir])
