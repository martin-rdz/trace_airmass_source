# probably docker ignore needs to be configured
docker build -t trace_env .

docker run --name trace_run -v /home/radenz/trace_pub/trace/..:/trace -it trace_env /bin/bash

docker image list

docker start trace_run
docker exec -it trace_run /bin/bash



docker run -it --rm -p 8890:8890 -v ~:/workingdir trace_env /bin/bash

jupyter notebook --ip=0.0.0.0 --no-browser --port=8890 --allow-root
