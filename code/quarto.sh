#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "Activate the conda environment first: conda activate differential-interaction-analysis" >&2
  exit 1
fi

export QUARTO_SHARE_PATH="${QUARTO_SHARE_PATH:-${CONDA_PREFIX}/share/quarto}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

export TMPDIR="${QUARTO_PROJECT_TMPDIR:-${PWD}/.quarto/tmp}"
mkdir -p "${TMPDIR}"

if [[ -z "${QUARTO_DART_SASS:-}" && -x "${CONDA_PREFIX}/bin/sass" ]]; then
  export QUARTO_DART_SASS="${CONDA_PREFIX}/bin/sass"
fi

arch="$(uname -m)"
case "${arch}" in
  arm64) deno_arch="aarch64" ;;
  x86_64) deno_arch="x86_64" ;;
  *) deno_arch="${arch}" ;;
esac

deno_dir="${CONDA_PREFIX}/bin/tools/${deno_arch}"
deno_link="${deno_dir}/deno"

if [[ ! -x "${deno_link}" && -x "${CONDA_PREFIX}/bin/deno" ]]; then
  mkdir -p "${deno_dir}"
  ln -sf ../../deno "${deno_link}"
fi

pandoc_link="${deno_dir}/pandoc"
if [[ ! -x "${pandoc_link}" && -x "${CONDA_PREFIX}/bin/pandoc" ]]; then
  mkdir -p "${deno_dir}"
  ln -sf ../../pandoc "${pandoc_link}"
fi

esbuild_link="${deno_dir}/esbuild"
if [[ ! -x "${esbuild_link}" && -x "${CONDA_PREFIX}/bin/esbuild" ]]; then
  mkdir -p "${deno_dir}"
  ln -sf ../../esbuild "${esbuild_link}"
fi

exec "${CONDA_PREFIX}/bin/quarto" "$@"
