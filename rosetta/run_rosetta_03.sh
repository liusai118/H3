#!/usr/bin/env bash
set -euo pipefail

# =========================
# =========================
ROSETTA_BIN="/home/liusai/software/rosetta/rosetta.binary.ubuntu.release-371/main/source/bin"
RELAX_EXE="${ROSETTA_BIN}/relax.static.linuxgccrelease"
IA_EXE="${ROSETTA_BIN}/InterfaceAnalyzer.static.linuxgccrelease"
BTK_PARAMS="/home/liusai/rosetta/ref/BTK/BTK.params"

ROOT_DIR="../structure"
OUT_ROOT="./rosetta_batch_out"
NSTRUCT=20
MAX_JOBS=8

SHOW_WARNINGS=1

# =========================
# =========================
ROOT_DIR="$(readlink -f "${ROOT_DIR}")"
OUT_ROOT="$(readlink -f "${OUT_ROOT}")"
RELAX_EXE="$(readlink -f "${RELAX_EXE}")"
IA_EXE="$(readlink -f "${IA_EXE}")"
BTK_PARAMS="$(readlink -f "${BTK_PARAMS}")"

mkdir -p "${OUT_ROOT}"

# =========================
# =========================
[[ -x "${RELAX_EXE}" ]] || { echo "ERROR: RELAX_EXE not found: ${RELAX_EXE}"; exit 1; }
[[ -x "${IA_EXE}"    ]] || { echo "ERROR: IA_EXE not found: ${IA_EXE}"; exit 1; }
[[ -f "${BTK_PARAMS}" ]] || { echo "ERROR: BTK_PARAMS not found: ${BTK_PARAMS}"; exit 1; }

log_warn() {
    if [[ "${SHOW_WARNINGS}" -eq 1 ]]; then
        echo "$*"
    fi
}

run_one_complex() {
    local complex_dir="$1"
    local complex_id
    complex_id=$(basename "${complex_dir}")

    local complex_out="${OUT_ROOT}/${complex_id}"
    mkdir -p "${complex_out}"

    mapfile -t pdb_files < <(
        find "${complex_dir}" -maxdepth 1 -type f -regextype posix-extended \
        -regex '.*/sample_[0-9]+\.pdb' | sort
    )

    if [[ ${#pdb_files[@]} -eq 0 ]]; then
        log_warn "[${complex_id}] No sample pdb files found, skip."
        return 0
    fi

    for pdb in "${pdb_files[@]}"; do
        pdb="$(readlink -f "${pdb}")"

        local base sample_out relax_dir interface_dir
        base=$(basename "${pdb}" .pdb)

        sample_out="${complex_out}/${base}"
        relax_dir="${sample_out}/relax"
        interface_dir="${sample_out}/interface"

        mkdir -p "${relax_dir}" "${interface_dir}"

        echo "[${complex_id}] Relax: ${base}"

        (
            cd "${relax_dir}"
            "${RELAX_EXE}" \
              -s "${pdb}" \
              -extra_res_fa "${BTK_PARAMS}" \
              -relax:fast \
              -relax:constrain_relax_to_start_coords \
              -relax:ramp_constraints false \
              -nstruct "${NSTRUCT}" \
              -out:path:all "${relax_dir}" \
              -out:file:scorefile "${relax_dir}/relax.sc" \
              -ex1 -ex2aro \
              -use_input_sc \
              -overwrite \
              > "${relax_dir}/relax.log" 2>&1
        ) || {
            log_warn "[${complex_id}] ERROR: relax failed for ${base} (see ${relax_dir}/relax.log)"
            continue
        }

        shopt -s nullglob
        local decoys=("${relax_dir}/${base}"_*.pdb)
        shopt -u nullglob

        if [[ ${#decoys[@]} -eq 0 ]]; then
            log_warn "[${complex_id}] WARNING: no decoys for ${base} (see ${relax_dir}/relax.log)"
            continue
        fi

        for decoy in "${decoys[@]}"; do
            decoy="$(readlink -f "${decoy}")"

            local decoy_base ia_log ia_score
            decoy_base=$(basename "${decoy}" .pdb)
            ia_log="${interface_dir}/${decoy_base}.interface.log"
            ia_score="${interface_dir}/${decoy_base}.interface.sc"

            echo "[${complex_id}] Interface: ${decoy_base}"

            (
                cd "${interface_dir}"
                "${IA_EXE}" \
                  -s "${decoy}" \
                  -extra_res_fa "${BTK_PARAMS}" \
                  -pack_input true \
                  -pack_separated true \
                  -add_regular_scores_to_scorefile true \
                  -tracer_data_print true \
                  -out:path:all "${interface_dir}" \
                  -out:file:score_only "${ia_score}" \
                  -overwrite \
                  > "${ia_log}" 2>&1
            ) || {
                log_warn "[${complex_id}] WARNING: interface failed for ${decoy_base} (see ${ia_log})"
                continue
            }

            if [[ -f "${ia_score}" ]]; then
                :
            elif grep -q "reported success" "${ia_log}" 2>/dev/null; then
                log_warn "[${complex_id}] WARNING: interface ran but no score file for ${decoy_base} (see ${ia_log})"
            else
                log_warn "[${complex_id}] WARNING: interface may be incomplete for ${decoy_base} (see ${ia_log})"
            fi
        done

        echo "[${complex_id}] Done: ${base}"
    done
}

export RELAX_EXE IA_EXE BTK_PARAMS OUT_ROOT NSTRUCT SHOW_WARNINGS
export -f run_one_complex
export -f log_warn

dirs=()
for d in "${ROOT_DIR}"/*; do
    [[ -d "$d" ]] || continue
    [[ "$(basename "$d")" == "$(basename "$OUT_ROOT")" ]] && continue
    [[ "$(basename "$d")" =~ ^[0-9]+$ ]] || continue
    dirs+=("$(readlink -f "$d")")
done

if [[ ${#dirs[@]} -eq 0 ]]; then
    echo "No complex directories found."
    exit 1
fi

echo "Found ${#dirs[@]} complex directories."
echo "MAX_JOBS=${MAX_JOBS}, NSTRUCT=${NSTRUCT}"

running=0
for d in "${dirs[@]}"; do
    run_one_complex "$d" &
    ((running+=1))

    if [[ "${running}" -ge "${MAX_JOBS}" ]]; then
        wait -n || true
        ((running-=1))
    fi
done

wait || true

echo "All relax + interface jobs finished."
echo "Output root: ${OUT_ROOT}"

