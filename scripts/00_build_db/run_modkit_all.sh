#!/bin/bash

# ---------- CONFIG ----------
SHEET="metadata/sample_sheet.tsv"
CH3_DIR="data/ch3"
OUT_DIR="results/bedmethyl"
LOG_DIR="logs/modkit"

/bin/mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

# ---------- LOOP ----------
tail -n +2 "$SHEET" | while IFS=$'\t' read -r barcode type notes; do
    CH3_FILE="${CH3_DIR}/${barcode}-0.ch3"
    OUT_FILE="${OUT_DIR}/${barcode}.bedmethyl"

    echo ">> Starting $barcode"

    if [ ! -f "$CH3_FILE" ]; then
        echo "   !! CH3 not found: $CH3_FILE"
        continue
    fi

    # launch modkit pileup (ruta ABSOLUTA)
/Users/antonioescudero/bin/modkit pileup \
        --only-cpg \
        --cpg-weighting none \
        "$CH3_FILE" \
        "$OUT_FILE" \
        > "${LOG_DIR}/${barcode}.log" 2>&1 &

    echo "   launched (log: ${LOG_DIR}/${barcode}.log)"
done

echo ">> All jobs launched in background."

