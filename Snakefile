rule meltingpoint:
    input:
        "input.json"
    output:
        "output.json",
        "script.nbconvert.ipynb",
        "plot.nbconvert.ipynb"
    conda:
        "envs/melting.yaml"
    shell:
        "export PYIRONRESOURCEPATHS='./resources'; export PYIRONPROJECTPATHS='.'; papermill ./scripts/script.ipynb script.nbconvert.ipynb -k 'python3' -p input_file 'input.json' -p output_file 'output.json'; papermill ./scripts/plot.ipynb plot.nbconvert.ipynb -k 'python3' -p input_file 'output.json'"
