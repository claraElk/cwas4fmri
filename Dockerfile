FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS build

# Install dependencies before the package itself to leverage caching
COPY . /app

ENV UV_NO_DEV=1
WORKDIR /app
RUN --mount=type=cache,target=/root/.cache/uv \
uv sync --locked --python 3.12

ENTRYPOINT ["uv", "run", "cwas4fmri"]
