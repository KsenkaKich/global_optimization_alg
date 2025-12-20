import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def testFunction1(x):
    return (x - 2.0) ** 2 + 1.0

def testFunction2(x):
    return x ** 2 + np.cos(18 * x)

def testFunction3(x):
    return np.exp(x) + np.sin(17 * x)

def visualize_algorithm(filename, algorithm_name, func):
    
    try:
        data = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"File {filename} not found")
        return
    
    x_min = data['x'].min() - 0.2
    x_max = data['x'].max() + 0.2
    
    x_vals = np.linspace(x_min, x_max, 500)
    y_vals = func(x_vals)
    
    #график функции
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_vals, y_vals, 'b-', linewidth=2, alpha=0.7)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('f(x)', fontsize=12)
    ax.set_title(f'{algorithm_name}', fontsize=14, fontweight='bold')
    
    y_min = y_vals.min() - 5
    y_max = y_vals.max() + 5
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    #основная сетка графика
    ax.grid(True, alpha=0.9, linestyle=':', linewidth=0.3)
    
    #количество делений на осях
    ax.xaxis.set_major_locator(plt.MaxNLocator(20)) 
    ax.yaxis.set_major_locator(plt.MaxNLocator(20))
    
    #вспомогательная мелкая сетка(минорная)
    ax.minorticks_on()
    ax.grid(which='minor', alpha=0.7, linestyle=':', linewidth=0.2)
    
    #лучшая точка
    best_idx = data['f(x)'].idxmin()
    best_x = data.loc[best_idx, 'x']
    best_y = data.loc[best_idx, 'f(x)']
    
    #все точки алгоритма
    ax.scatter(data['x'], data['f(x)'], color='purple', s=30, 
               edgecolor='black', linewidth=0.5, zorder=5, alpha=0.8)
    
    #лучшая точка - звёздочка
    ax.scatter(best_x, best_y, color='red', s=100, marker='*', 
               edgecolor='black', linewidth=1, zorder=10)
        
    #избегание наложений
    plt.tight_layout()
    
    #статистика
    print(f"\n{algorithm_name}:")
    print(f"Итераций: {len(data)}")
    print(f"Лучшая точка: x = {best_x:.6f}, f(x) = {best_y:.6f}")
    
    output_filename = filename.replace('.csv', '_graphic.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    
    # временная мера
    current_func = testFunction3
    
    visualize_algorithm("gsa_results.csv", "GSA Algoritm", current_func)
    visualize_algorithm("scan_results.csv", "Scan Algoritm", current_func)