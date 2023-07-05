# Fragmenstein workflows

This repository is used to build the things needed to execute
the **Squonk2 Fragmenstein** jobs.

## Building the workflows
GitHib Actions take care of official releases of these Job workflows. A container
image is built and pushed to DockerHub using the image tag `latest` for each change
on the main branch. Container images tagged with an official release (e.g. using the
standard format `1.0.0`) will be built and pushed using the tag used for the release.

Refer to the `docker-compose.yaml` file for the details of the images that are built.

If you want to build your own images you can use the docker compose file in the
repository, but **DO NOT** publish your own images to DockerHub: -

    python -m venv venv
    source venv/bin/activate
    pip install -r build-requirements.txt

    docker-compose build
