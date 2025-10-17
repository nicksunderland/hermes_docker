docker run --rm \
  -v /Users/xx20081/git/hermes_docker/tests/data:/data \
  -v /Users/xx20081/git/hermes_docker/tests/resources:/resources \
  -v /Users/xx20081/git/hermes_docker/tests/output:/output \
  nicksunderland/hermes_docker:latest \
  /data/test_data.tsv.gz \
  /data/config.json \
  /resources \
  /output


# server command ???
#ContainerProperties:
#  Image: nicksunderland/hermes_docker:latest
#  Command:
#    - '/bin/bash'
#    - '/opt/run.sh'
#    - 'Ref::s3-path'           # $1 → GWAS input (S3 path)
#    - 'Ref::config-json'       # $2 → JSON file or inline JSON string (the gwas meta-data)
#    - 'Ref::s3-resources-dir'  # $3 → reference resources dir (S3 hermes repository directory)
#    - 'Ref::file-guid'         # $4 → GUID of file
