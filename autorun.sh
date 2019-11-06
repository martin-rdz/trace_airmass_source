#!/bin/bash

#python3 gen_hysplit_input.py --station punta
daterange=20181221-20190121
daterange=20190121-20190201
#daterange=20190201-20190203
python3 test_data_avail.py --station punta --daterange $daterange
python3 run_assemble_hysplit.py --station punta --daterange $daterange
python3 plot2d.py --station punta --daterange $daterange
#python3 compress_data.py --station punta --daterange 20181210-20181221

#python3 test_data_avail.py --station ps106 --daterange 20170525-20170719
#python3 run_assemble_hysplit.py --station ps106 --daterange 20170525-20170719
#python3 run_assemble_hysplit.py --station ps106 --daterange 20170620-20170723
#python3 plot2d.py --station ps106 --daterange 20170525-20170719
#python3 compress_data.py --station ps116 --daterange 20181111-20181211
