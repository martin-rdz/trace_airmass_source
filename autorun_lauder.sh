#!/bin/bash

set -e

#python3 gen_hysplit_input.py --station punta

daterange=20180101-20180111
daterange=20180111-20180131

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange

#daterange=20180201-20180211
daterange=20180211-20180228

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange

daterange=20180611-20180630

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange


daterange=20180301-20180311
daterange=20180311-20180331

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange

daterange=20180401-20180411
daterange=20180411-20180430

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange

daterange=20180701-20180711
daterange=20180711-20180731

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange

daterange=20180501-20180531

#python3 test_data_avail.py --station lauder --daterange $daterange
#python3 run_assemble_hysplit.py --station lauder --daterange $daterange
#python3 plot2d.py --station lauder --daterange $daterange
#python3 compress_data.py --station lauder --daterange $daterange

python3 plot2d.py --station lauder --daterange 20180112-20180731
