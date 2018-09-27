#! /usr/bin/env python3
# coding=utf-8
"""
Author: radenz@tropos.de

"""

import datetime

import trace_source.hysplit
import trace_source.land_sfc
import trace_source.helpers
import trace_source.simulations


def time_list(begin, end, delta):
    """generate a time list from begin to <= end with delta"""
    time_list = []
    elem = begin
    while elem <= end:
        time_list.append(elem)
        elem += datetime.timedelta(hours=delta)
    return time_list

