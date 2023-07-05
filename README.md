# Fragmenstein workflows

![GitHub](https://img.shields.io/github/license/InformaticsMatters/squonk2-fragmenstein)

![GitHub release (with filter)](https://img.shields.io/github/v/release/InformaticsMatters/squonk2-fragmenstein)

[![latest](https://github.com/InformaticsMatters/squonk2-fragmenstein/actions/workflows/publish-latest.yaml/badge.svg)](https://github.com/InformaticsMatters/squonk2-fragmenstein/actions/workflows/publish-latest.yaml)
[![stable](https://github.com/InformaticsMatters/squonk2-fragmenstein/actions/workflows/publish-stable.yaml/badge.svg)](https://github.com/InformaticsMatters/squonk2-fragmenstein/actions/workflows/publish-stable.yaml)
[![build](https://github.com/InformaticsMatters/squonk2-fragmenstein/actions/workflows/build.yaml/badge.svg)](https://github.com/InformaticsMatters/squonk2-fragmenstein/actions/workflows/build.yaml)

This repository is used to build the things needed to execute
the **Squonk2 Fragmenstein** jobs.

## Building the workflows
Automated GitHub CI Actions take care of official releases of these Job workflows.

Container images are built and pushed to DockerHub using the image tag `latest` for
each change on the main branch, and images built during an _official release_
will be tagged `stable`.

>   To build an _official release_ create a new **Release** in GitHub,
    assigning a [semver] style tag to the release, while also leaving the
    **Set as a pre-release** checkbox un-ticked.

Refer to the `docker-compose.yaml` file for the details of the images that are built.

If you want to build your own (local) images you can use the docker compose file in the
repository, but **DO NOT** publish these images to DockerHub: -

    python -m venv venv
    source venv/bin/activate
    pip install -r build-requirements.txt

    docker-compose build

---

[semver]: https://semver.org
