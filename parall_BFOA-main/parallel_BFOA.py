import time
from copy import deepcopy
import matplotlib.pyplot as plt
from bacteria import bacteria
from fastaReader import fastaReader
import random

# SIN Manager, listas normales para aligerar

def run_bfoa(
    dAttr=0.05, wAttr=0.005, hRep=None, wRep=0.01,
    iteraciones=30, numeroDeBacterias=10, tumbo=2
):
    if hRep is None:
        hRep = dAttr

    secuencias = fastaReader().seqs
    names = fastaReader().names

    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])

    numSec = len(secuencias)

    # USAR listas normales, NO Manager
    poblacion = [None for _ in range(numeroDeBacterias)]
    names = list(names)

    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)

    def mutacion_adaptativa(poblacion, bestIdx, intensidad=1):
        for i in range(len(poblacion)):
            if i != bestIdx:
                bacterium = list(poblacion[i])
                for _ in range(intensidad):
                    seq_idx = random.randint(0, len(bacterium) - 1)
                    pos = random.randint(0, len(bacterium[seq_idx]))
                    bacterium[seq_idx] = bacterium[seq_idx][:pos] + ["-"] + bacterium[seq_idx][pos:]
                poblacion[i] = bacterium  # Dejarlas como lista

    operadorBacterial = bacteria(numeroDeBacterias)
    veryBest = [None, None, None]
    globalNFE = 0
    fitness_history = []

    start_time = time.time()
    poblacionInicial()

    print(f"\n=== Ejecutando prueba con dAttr={dAttr}, wAttr={wAttr}, hRep={hRep}, wRep={wRep}, iteraciones={iteraciones}, bacterias={numeroDeBacterias}, tumbo={tumbo} ===")

    for it in range(iteraciones):
        print(f"\n-- Iteraci\u00f3n {it + 1} --")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)
        operadorBacterial.cuadra(numSec, poblacion)
        operadorBacterial.creaGranListaPares(poblacion)
        operadorBacterial.evaluaBlosum()
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        operadorBacterial.creaTablaInteraction()
        operadorBacterial.creaTablaFitness()

        globalNFE += operadorBacterial.getNFE()
        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
        fitness_history.append(bestFitness)

        print(f"Best ID: {bestIdx} | Fitness: {operadorBacterial.tablaFitness[bestIdx]:.2f} | "
              f"Blosum: {operadorBacterial.blosumScore[bestIdx]:.2f} | "
              f"Interaction: {operadorBacterial.tablaInteraction[bestIdx]:.2f} | NFE total: {globalNFE}")

        if (veryBest[0] is None) or (bestFitness > veryBest[1]):
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = deepcopy(poblacion[bestIdx])

        mutacion_adaptativa(poblacion, bestIdx, intensidad=1)

        operadorBacterial.replaceWorst(poblacion, veryBest[0])
        operadorBacterial.resetListas(numeroDeBacterias)

    total_time = time.time() - start_time
    print(f"\n✅ Mejor fitness global: {veryBest[1]:.2f}")
    print(f"⏱️ Tiempo total: {total_time:.2f} segundos\n")

    # 📈 Gr\u00e1fica de evoluci\u00f3n del fitness
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, iteraciones + 1), fitness_history, marker='o', color='blue', label='Fitness por iteraci\u00f3n')
    plt.title("Evoluci\u00f3n del Fitness por Iteraci\u00f3n")
    plt.xlabel("Iteraci\u00f3n")
    plt.ylabel("Fitness")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return veryBest[1], total_time


if __name__ == "__main__":
    try:
        print("\n=== Ejecutando en modo ligero ===")
        run_bfoa()
    except Exception as e:
        print("❌ Error durante la ejecuci\u00f3n:", str(e))
