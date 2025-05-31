import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

images_path = "images/"


# Função para extrair dados da tabela final de um ficheiro de texto
def extrair_dados_tune(content):
    # Regex com grupos opcionais a partir do "Tempo Médio (s)"
    pattern = (
        r"(\d+)\s*\|\s*"  # r
        r"([\d.]+)\s*\|\s*"  # Min SP
        r"([\d.]+)\s*\|\s*"  # Média SP
        r"([\d.]+)"  # Max SP
        r"(?:\s*\|\s*([\d.]+))?"  # Tempo Médio (s) (opcional)
        r"(?:\s*\|\s*([\d.]+))?"  # Iterações GRASP (opcional)
        r"(?:\s*\|\s*([\d.]+))?"  # Iterações Local Search (opcional)
    )

    matches = re.findall(pattern, content)
    df = pd.DataFrame(
        matches,
        columns=[
            "r",
            "Min SP",
            "Média SP",
            "Max SP",
            "Tempo Médio (s)",
            "Iterações GRASP",
            "Iterações Local Search",
        ],
    )

    # Converte para float e int, tratando campos ausentes como NaN
    df = df.replace("", np.nan)
    df = df.astype(
        {
            "r": int,
            "Min SP": float,
            "Média SP": float,
            "Max SP": float,
            "Tempo Médio (s)": float,
            "Iterações GRASP": float,
            "Iterações Local Search": float,
        }
    )

    return df


# Leitura dos ficheiros (substituir os nomes se estiverem em caminhos diferentes)
with open("tune_10runs_30sec.txt", "r") as f1:
    content_10_30 = f1.read()

with open("tune_30runs_10sec.txt", "r") as f2:
    content_30_10 = f2.read()

with open("tune_10runs_30sec_optimized.txt", "r") as f3:
    content_10_30_opt = f3.read()

with open("tune_30runs_10sec_optimized.txt", "r") as f4:
    content_30_10_opt = f4.read()

# Extrair dataframes
df_10_30 = extrair_dados_tune(content_10_30)
df_30_10 = extrair_dados_tune(content_30_10)
df_10_30_opt = extrair_dados_tune(content_10_30_opt)
df_30_10_opt = extrair_dados_tune(content_30_10_opt)

# Comparação Otimizada vs Não Otimizada para 30 runs e 10 segundos
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(df_10_30["r"], df_10_30["Média SP"], marker="o", label="Não Otimizada")
axes[0].plot(df_10_30_opt["r"], df_10_30_opt["Média SP"], marker="o", label="Otimizada")
axes[0].set_title("Média do Shortest Path (SP) vs r")
axes[0].set_xlabel("Parâmetro r")
axes[0].set_ylabel("Média SP")
axes[0].grid(True)
axes[0].legend()

axes[1].plot(
    df_10_30["r"], df_10_30["Tempo Médio (s)"], marker="o", label="Não Otimizada"
)
axes[1].plot(
    df_10_30_opt["r"], df_10_30_opt["Tempo Médio (s)"], marker="o", label="Otimizada"
)
axes[1].set_title("Tempo Médio de Execução vs r")
axes[1].set_xlabel("Parâmetro r")
axes[1].set_ylabel("Tempo Médio (s)")
axes[1].grid(True)
axes[1].legend()

plt.tight_layout()
plt.savefig(images_path + "comparacao_sp_tempo_opt_vs_normal.png")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(
    df_10_30["r"], df_10_30["Iterações GRASP"], marker="o", label="Não Otimizada"
)
axes[0].plot(
    df_10_30_opt["r"], df_10_30_opt["Iterações GRASP"], marker="o", label="Otimizada"
)
axes[0].set_title("Iterações GRASP vs r")
axes[0].set_xlabel("Parâmetro r")
axes[0].set_ylabel("Iterações GRASP")
axes[0].grid(True)
axes[0].legend()

axes[1].plot(
    df_10_30["r"], df_10_30["Iterações Local Search"], marker="o", label="Não Otimizada"
)
axes[1].plot(
    df_10_30_opt["r"],
    df_10_30_opt["Iterações Local Search"],
    marker="o",
    label="Otimizada",
)
axes[1].set_title("Iterações Local Search vs r")
axes[1].set_xlabel("Parâmetro r")
axes[1].set_ylabel("Iterações Local Search")
axes[1].grid(True)
axes[1].legend()

plt.tight_layout()
plt.savefig(images_path + "comparacao_iterations_opt_vs_normal.png")

# Comparação entre 10 runs de 30 segundos e 30 runs de 10 segundos (sem otimização)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(df_10_30["r"], df_10_30["Média SP"], marker="o", label="10s / 30 runs")
axes[0].plot(df_30_10["r"], df_30_10["Média SP"], marker="o", label="30s / 10 runs")
axes[0].set_title("Média do Shortest Path (SP) vs r")
axes[0].set_xlabel("Parâmetro r")
axes[0].set_ylabel("Média SP")
axes[0].grid(True)
axes[0].legend()

axes[1].plot(
    df_10_30["r"], df_10_30["Tempo Médio (s)"], marker="o", label="10s / 30 runs"
)
axes[1].plot(
    df_30_10["r"], df_30_10["Tempo Médio (s)"], marker="o", label="30s / 10 runs"
)
axes[1].set_title("Tempo Médio de Execução vs r")
axes[1].set_xlabel("Parâmetro r")
axes[1].set_ylabel("Tempo Médio (s)")
axes[1].grid(True)
axes[1].legend()

plt.tight_layout()
plt.savefig(images_path + "comparacao_10_30_vs_30_10_runs.png")

# Comparação entre 10 runs de 30 segundos e 30 runs de 10 segundos (com otimização)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].plot(
    df_10_30_opt["r"], df_10_30_opt["Média SP"], marker="o", label="10s / 30 runs"
)
axes[0].plot(
    df_30_10_opt["r"], df_30_10_opt["Média SP"], marker="o", label="30s / 10 runs"
)
axes[0].set_title("Média do Shortest Path (SP) vs r (Otimizada)")
axes[0].set_xlabel("Parâmetro r")
axes[0].set_ylabel("Média SP")
axes[0].grid(True)
axes[0].legend()

axes[1].plot(
    df_10_30_opt["r"],
    df_10_30_opt["Tempo Médio (s)"],
    marker="o",
    label="10s / 30 runs",
)
axes[1].plot(
    df_30_10_opt["r"],
    df_30_10_opt["Tempo Médio (s)"],
    marker="o",
    label="30s / 10 runs",
)
axes[1].set_title("Tempo Médio de Execução vs r (Otimizada)")
axes[1].set_xlabel("Parâmetro r")
axes[1].set_ylabel("Tempo Médio (s)")
axes[1].grid(True)
axes[1].legend()

plt.tight_layout()
plt.savefig(images_path + "comparacao_10_30_vs_30_10_runs_opt.png")
