# Cassandra Bodin
# Phys305
# Final Project: Kirkwood Gap Explorer
# UI tool for exploring parameter space and data mining of asteroid belt simulations

import tkinter as tk
from tkinter import ttk, messagebox
import sqlite3
import multiprocessing
import time
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
import astropy.units as u

# Global constants
G = 6.674e-11
GMs = 4 * (math.pi ** 2)
m_sun = 1.989e30

def create_database():
    """Create SQLite database to store simulation results"""
    conn = sqlite3.connect('kirkwood_gap_simulations.db')
    cursor = conn.cursor()
    
    # Create simulations table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS simulations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            jupiter_mass REAL,
            jupiter_ecc REAL,
            num_asteroids INTEGER,
            simulation_time REAL,
            start_time TEXT,
            end_time TEXT,
            status TEXT
        )
    ''')
    
    # Create results table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            simulation_id INTEGER,
            semi_major_axes TEXT,
            FOREIGN KEY (simulation_id) REFERENCES simulations(id)
        )
    ''')
    
    conn.commit()
    conn.close()

def run_simulation(params):
    """Run a single simulation with given parameters"""
    jupiter_mass, jupiter_ecc, num_asteroids, simulation_time = params
    
    try:
        # Initialize Jupiter with given parameters
        solar_system_ephemeris.set('jpl')
        t = Time('J2000', format='jyear')
        jupiter_pos, jupiter_vel = get_body_barycentric_posvel('jupiter', t)
        
        # Convert to AU and AU/year
        x2 = jupiter_pos.x.to(u.au).value
        y2 = jupiter_pos.y.to(u.au).value
        vx2 = jupiter_vel.x.to(u.au/u.year).value
        vy2 = jupiter_vel.y.to(u.au/u.year).value
        
        # Generate asteroids
        n = 1.5 * np.random.random(num_asteroids) + 2.0  # Semi-major axes between 2 and 3.5 AU
        m_a1 = 1e10  # Mass of asteroids (arbitrary)
        
        x = np.array([])
        y = np.array([])
        vx = np.array([])
        vy = np.array([])
        semimajor_a = np.array([])
        
        for semimajor_a1 in n:
            a_a1 = semimajor_a1
            ecc_a1 = 0  # Circular orbits for asteroids
            perihelion_a1 = a_a1 * (1 - ecc_a1)
            
            # Calculate initial velocity for circular orbit
            vp_a1 = math.sqrt(GMs) * math.sqrt((1 + ecc_a1) / (a_a1 * (1 - ecc_a1)) * (1 + (m_a1 / m_sun)))
            astr_vx = 0
            astr_vy = vp_a1
            astr_x = -perihelion_a1
            astr_y = 0
            
            x = np.append(x, astr_x)
            y = np.append(y, astr_y)
            vx = np.append(vx, astr_vx)
            vy = np.append(vy, astr_vy)
            semimajor_a = np.append(semimajor_a, semimajor_a1)
        
        # Simulation parameters
        h = 1e-2  # Time step in years
        t_initial = 0
        t = t_initial
        count = 0
        
        # Variables for storing positions
        xarr = np.array([])
        yarr = np.array([])
        xarr = np.append(xarr, x)
        yarr = np.append(yarr, y)
        
        xarr2 = np.array([])
        yarr2 = np.array([])
        xarr2 = np.append(xarr2, x2)
        yarr2 = np.append(yarr2, y2)
        
        # Run simulation
        while t < simulation_time:
            # Calculate radial distances
            r2 = math.sqrt(x2 ** 2 + y2 ** 2)
            
            # Update asteroid velocities
            for i in range(num_asteroids):
                r = math.sqrt(x[i] ** 2 + y[i] ** 2)
                rrel = math.sqrt((x[i] - x2) ** 2 + (y[i] - y2) ** 2)
                
                vx[i] += h * (-GMs * x[i] / r ** 3 - GMs * (jupiter_mass / m_sun) * (x[i] - x2) / rrel ** 3)
                vy[i] += h * (-GMs * y[i] / r ** 3 - GMs * (jupiter_mass / m_sun) * (y[i] - y2) / rrel ** 3)
            
            # Update Jupiter velocity
            vx2 += h * (-GMs * x2 / r2 ** 3)
            vy2 += h * (-GMs * y2 / r2 ** 3)
            
            # Update asteroid positions
            for i in range(num_asteroids):
                x[i] += vx[i] * h
                y[i] += vy[i] * h
            
            # Update Jupiter position
            x2 += vx2 * h
            y2 += vy2 * h
            
            # Store positions
            xarr = np.append(xarr, x)
            yarr = np.append(yarr, y)
            xarr2 = np.append(xarr2, x2)
            yarr2 = np.append(yarr2, y2)
            
            count += 1
            t += h
        
        # Reshape position arrays
        xarr = xarr.reshape([num_asteroids, count + 1])
        yarr = yarr.reshape([num_asteroids, count + 1])
        
        # Calculate final semi-major axes
        r_final = np.sqrt((xarr[:, -1]) ** 2 + (yarr[:, -1]) ** 2)
        
        # Filter asteroids that are still in the belt (2-3.5 AU)
        surviving_asteroids = semimajor_a[(r_final >= 2.0) & (r_final <= 3.5)]
        
        return surviving_asteroids
        
    except Exception as e:
        print(f"Simulation failed: {e}")
        return None

def calculate_gap_clarity(semi_major_axes):
    """Calculate the clarity of Kirkwood gaps"""
    if len(semi_major_axes) == 0:
        return 0.0
    
    # Define resonance locations (in AU) for 3:1, 5:2, 7:3, and 2:1 resonances
    resonances = [2.5, 2.8, 3.0, 3.3]
    
    # Create histogram
    bins = np.linspace(2.0, 3.5, 100)
    hist, bin_edges = np.histogram(semi_major_axes, bins=bins)
    
    # Calculate gap clarity for each resonance
    gap_clarities = []
    for res in resonances:
        # Find the bin containing the resonance
        bin_idx = np.digitize(res, bin_edges) - 1
        
        # Skip if resonance is outside the bin range
        if bin_idx < 0 or bin_idx >= len(hist):
            continue
        
        # Calculate the average density in neighboring bins
        window_size = 5
        start = max(0, bin_idx - window_size)
        end = min(len(hist), bin_idx + window_size + 1)
        
        # Exclude the resonance bin itself
        neighbor_bins = np.concatenate((hist[start:bin_idx], hist[bin_idx+1:end]))
        
        if len(neighbor_bins) == 0:
            continue
            
        avg_neighbor_density = np.mean(neighbor_bins)
        resonance_density = hist[bin_idx]
        
        # Calculate gap clarity as the ratio of resonance density to average neighbor density
        # Lower values indicate clearer gaps
        if avg_neighbor_density > 0:
            gap_clarity = resonance_density / avg_neighbor_density
            gap_clarities.append(gap_clarity)
    
    # Return the average gap clarity across all resonances
    if len(gap_clarities) == 0:
        return 1.0  # No gaps detected, maximum clarity value
    
    return np.mean(gap_clarities)

def calculate_survival_rate(semi_major_axes, initial_count):
    """Calculate the survival rate of asteroids"""
    if initial_count == 0:
        return 0.0
    
    return (len(semi_major_axes) / initial_count) * 100

def batch_simulator_tab(parent):
    """Create the batch simulator tab"""
    frame = ttk.Frame(parent)
    
    # Parameter grid
    ttk.Label(frame, text="木星质量范围 (太阳质量倍数)").grid(row=0, column=0, padx=5, pady=5, sticky="w")
    jup_mass_min = ttk.Entry(frame, width=10)
    jup_mass_min.grid(row=0, column=1, padx=5, pady=5)
    jup_mass_min.insert(0, "0.5")
    
    ttk.Label(frame, text="到").grid(row=0, column=2, padx=5, pady=5)
    
    jup_mass_max = ttk.Entry(frame, width=10)
    jup_mass_max.grid(row=0, column=3, padx=5, pady=5)
    jup_mass_max.insert(0, "1.5")
    
    ttk.Label(frame, text="步长").grid(row=0, column=4, padx=5, pady=5)
    jup_mass_steps = ttk.Entry(frame, width=10)
    jup_mass_steps.grid(row=0, column=5, padx=5, pady=5)
    jup_mass_steps.insert(0, "10")
    
    ttk.Label(frame, text="木星偏心率范围").grid(row=1, column=0, padx=5, pady=5, sticky="w")
    jup_ecc_min = ttk.Entry(frame, width=10)
    jup_ecc_min.grid(row=1, column=1, padx=5, pady=5)
    jup_ecc_min.insert(0, "0.0")
    
    ttk.Label(frame, text="到").grid(row=1, column=2, padx=5, pady=5)
    
    jup_ecc_max = ttk.Entry(frame, width=10)
    jup_ecc_max.grid(row=1, column=3, padx=5, pady=5)
    jup_ecc_max.insert(0, "0.2")
    
    ttk.Label(frame, text="步长").grid(row=1, column=4, padx=5, pady=5)
    jup_ecc_steps = ttk.Entry(frame, width=10)
    jup_ecc_steps.grid(row=1, column=5, padx=5, pady=5)
    jup_ecc_steps.insert(0, "10")
    
    ttk.Label(frame, text="小行星数量").grid(row=2, column=0, padx=5, pady=5, sticky="w")
    num_asteroids = ttk.Entry(frame, width=10)
    num_asteroids.grid(row=2, column=1, padx=5, pady=5)
    num_asteroids.insert(0, "100")
    
    ttk.Label(frame, text="模拟时间 (年)").grid(row=2, column=2, padx=5, pady=5, sticky="w")
    sim_time = ttk.Entry(frame, width=10)
    sim_time.grid(row=2, column=3, padx=5, pady=5)
    sim_time.insert(0, "200")
    
    # Start button
    start_button = ttk.Button(frame, text="开始批处理模拟")
    start_button.grid(row=3, column=0, columnspan=6, pady=10)
    
    # Task queue window
    queue_frame = ttk.LabelFrame(frame, text="任务队列")
    queue_frame.grid(row=4, column=0, columnspan=6, padx=5, pady=5, sticky="nsew")
    
    queue_text = tk.Text(queue_frame, height=15, width=80)
    queue_text.grid(row=0, column=0, sticky="nsew")
    
    scrollbar = ttk.Scrollbar(queue_frame, orient="vertical", command=queue_text.yview)
    scrollbar.grid(row=0, column=1, sticky="ns")
    queue_text.configure(yscrollcommand=scrollbar.set)
    
    # Configure grid weights
    frame.grid_columnconfigure(5, weight=1)
    frame.grid_rowconfigure(4, weight=1)
    queue_frame.grid_columnconfigure(0, weight=1)
    queue_frame.grid_rowconfigure(0, weight=1)
    
    def start_batch_simulation():
        """Start the batch simulation process"""
        try:
            # Get parameters from UI
            mass_min = float(jup_mass_min.get())
            mass_max = float(jup_mass_max.get())
            mass_steps = int(jup_mass_steps.get())
            
            ecc_min = float(jup_ecc_min.get())
            ecc_max = float(jup_ecc_max.get())
            ecc_steps = int(jup_ecc_steps.get())
            
            num_ast = int(num_asteroids.get())
            sim_t = float(sim_time.get())
            
            # Validate parameters
            if mass_min < 0 or mass_max < mass_min or mass_steps < 1:
                raise ValueError("无效的木星质量参数")
            
            if ecc_min < 0 or ecc_max > 1 or ecc_max < ecc_min or ecc_steps < 1:
                raise ValueError("无效的木星偏心率参数")
            
            if num_ast < 1 or sim_t < 0:
                raise ValueError("无效的小行星数量或模拟时间")
            
            # Generate parameter grid
            mass_values = np.linspace(mass_min, mass_max, mass_steps)
            ecc_values = np.linspace(ecc_min, ecc_max, ecc_steps)
            
            total_simulations = mass_steps * ecc_steps
            
            queue_text.insert(tk.END, f"开始批处理模拟: {total_simulations} 次模拟\n")
            queue_text.insert(tk.END, f"木星质量范围: {mass_min} 到 {mass_max} ({mass_steps} 步)\n")
            queue_text.insert(tk.END, f"木星偏心率范围: {ecc_min} 到 {ecc_max} ({ecc_steps} 步)\n")
            queue_text.insert(tk.END, f"小行星数量: {num_ast}\n")
            queue_text.insert(tk.END, f"模拟时间: {sim_t} 年\n")
            queue_text.insert(tk.END, "=" * 60 + "\n")
            queue_text.see(tk.END)
            
            # Create database connection
            conn = sqlite3.connect('kirkwood_gap_simulations.db')
            cursor = conn.cursor()
            
            # Create list of simulation parameters
            params_list = []
            for mass in mass_values:
                for ecc in ecc_values:
                    # Insert simulation into database
                    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
                    cursor.execute('''
                        INSERT INTO simulations (jupiter_mass, jupiter_ecc, num_asteroids, simulation_time, start_time, status)
                        VALUES (?, ?, ?, ?, ?, ?)
                    ''', (mass, ecc, num_ast, sim_t, start_time, 'running'))
                    
                    simulation_id = cursor.lastrowid
                    params_list.append((simulation_id, mass, ecc, num_ast, sim_t))
            
            conn.commit()
            conn.close()
            
            # Run simulations in parallel
            def run_simulation_wrapper(params):
                simulation_id, mass, ecc, num_ast, sim_t = params
                result = run_simulation((mass, ecc, num_ast, sim_t))
                
                # Update database with results
                conn = sqlite3.connect('kirkwood_gap_simulations.db')
                cursor = conn.cursor()
                
                if result is not None:
                    # Convert array to string for storage
                    result_str = ','.join(map(str, result))
                    cursor.execute('''
                        INSERT INTO results (simulation_id, semi_major_axes)
                        VALUES (?, ?)
                    ''', (simulation_id, result_str))
                    
                    end_time = time.strftime("%Y-%m-%d %H:%M:%S")
                    cursor.execute('''
                        UPDATE simulations SET status = ?, end_time = ? WHERE id = ?
                    ''', ('completed', end_time, simulation_id))
                else:
                    end_time = time.strftime("%Y-%m-%d %H:%M:%S")
                    cursor.execute('''
                        UPDATE simulations SET status = ?, end_time = ? WHERE id = ?
                    ''', ('failed', end_time, simulation_id))
                
                conn.commit()
                conn.close()
                
                return simulation_id, result is not None
            
            # Use multiprocessing to run simulations
            with multiprocessing.Pool() as pool:
                results = pool.map(run_simulation_wrapper, params_list)
            
            # Update task queue
            completed = sum(1 for _, success in results if success)
            failed = sum(1 for _, success in results if not success)
            
            queue_text.insert(tk.END, f"\n批处理模拟完成!\n")
            queue_text.insert(tk.END, f"成功: {completed} 次\n")
            queue_text.insert(tk.END, f"失败: {failed} 次\n")
            queue_text.insert(tk.END, "=" * 60 + "\n")
            queue_text.see(tk.END)
            
            messagebox.showinfo("批处理完成", f"模拟已完成! 成功: {completed}, 失败: {failed}")
            
        except ValueError as e:
            messagebox.showerror("参数错误", str(e))
        except Exception as e:
            messagebox.showerror("错误", f"批处理模拟失败: {str(e)}")
    
    start_button.config(command=start_batch_simulation)
    
    return frame

def data_mining_tab(parent):
    """Create the data mining and visualization tab"""
    frame = ttk.Frame(parent)
    
    # Heatmap frame
    heatmap_frame = ttk.LabelFrame(frame, text="热力图")
    heatmap_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
    
    # Create figure and canvas for heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    canvas = FigureCanvasTkAgg(fig, master=heatmap_frame)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=0, column=0, sticky="nsew")
    
    # Control frame
    control_frame = ttk.LabelFrame(frame, text="控制面板")
    control_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")
    
    # Z-axis selection
    ttk.Label(control_frame, text="选择热力图Z轴:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
    z_axis_var = tk.StringVar()
    z_axis_combo = ttk.Combobox(control_frame, textvariable=z_axis_var, values=["柯克伍德带隙清晰度", "小行星平均存活率"])
    z_axis_combo.grid(row=0, column=1, padx=5, pady=5)
    z_axis_combo.current(0)
    
    # Refresh button
    refresh_button = ttk.Button(control_frame, text="刷新热力图")
    refresh_button.grid(row=1, column=0, columnspan=2, padx=5, pady=10)
    
    # Statistics frame
    stats_frame = ttk.LabelFrame(frame, text="统计信息")
    stats_frame.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="nsew")
    
    stats_text = tk.Text(stats_frame, height=10, width=100)
    stats_text.grid(row=0, column=0, sticky="nsew")
    
    # Configure grid weights
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    heatmap_frame.grid_columnconfigure(0, weight=1)
    heatmap_frame.grid_rowconfigure(0, weight=1)
    control_frame.grid_columnconfigure(1, weight=1)
    stats_frame.grid_columnconfigure(0, weight=1)
    stats_frame.grid_rowconfigure(0, weight=1)
    
    def update_heatmap():
        """Update the heatmap with current data"""
        try:
            conn = sqlite3.connect('kirkwood_gap_simulations.db')
            cursor = conn.cursor()
            
            # Get all completed simulations
            cursor.execute('''
                SELECT s.id, s.jupiter_mass, s.jupiter_ecc, s.num_asteroids, r.semi_major_axes
                FROM simulations s
                JOIN results r ON s.id = r.simulation_id
                WHERE s.status = 'completed'
            ''')
            simulations = cursor.fetchall()
            
            if not simulations:
                messagebox.showinfo("无数据", "数据库中没有已完成的模拟结果")
                return
            
            # Extract parameters and results
            jupiter_masses = []
            jupiter_eccs = []
            gap_clarities = []
            survival_rates = []
            
            for sim_id, mass, ecc, num_ast, semi_majors_str in simulations:
                # Convert semi-major axes string back to array
                semi_majors = np.array(list(map(float, semi_majors_str.split(','))))
                
                # Calculate gap clarity
                clarity = calculate_gap_clarity(semi_majors)
                
                # Calculate survival rate
                survival_rate = calculate_survival_rate(semi_majors, num_ast)
                
                jupiter_masses.append(mass)
                jupiter_eccs.append(ecc)
                gap_clarities.append(clarity)
                survival_rates.append(survival_rate)
            
            # Convert to numpy arrays
            jupiter_masses = np.array(jupiter_masses)
            jupiter_eccs = np.array(jupiter_eccs)
            gap_clarities = np.array(gap_clarities)
            survival_rates = np.array(survival_rates)
            
            # Create grid for heatmap
            mass_unique = np.unique(jupiter_masses)
            ecc_unique = np.unique(jupiter_eccs)
            
            # Create 2D arrays for heatmap data
            if z_axis_var.get() == "柯克伍德带隙清晰度":
                z_data = np.zeros((len(ecc_unique), len(mass_unique)))
                for i, ecc in enumerate(ecc_unique):
                    for j, mass in enumerate(mass_unique):
                        mask = (jupiter_masses == mass) & (jupiter_eccs == ecc)
                        if np.any(mask):
                            z_data[i, j] = gap_clarities[mask][0]
                
                cmap = 'viridis_r'  # Reverse viridis for clearer gaps (lower values are better)
                z_label = "柯克伍德带隙清晰度 (值越低越清晰)"
                title = "木星质量和偏心率对柯克伍德带隙的影响"
            else:
                z_data = np.zeros((len(ecc_unique), len(mass_unique)))
                for i, ecc in enumerate(ecc_unique):
                    for j, mass in enumerate(mass_unique):
                        mask = (jupiter_masses == mass) & (jupiter_eccs == ecc)
                        if np.any(mask):
                            z_data[i, j] = survival_rates[mask][0]
                
                cmap = 'viridis'
                z_label = "小行星平均存活率 (%)"
                title = "木星质量和偏心率对小行星存活率的影响"
            
            # Clear previous plot
            ax.clear()
            
            # Create heatmap
            im = ax.imshow(z_data, cmap=cmap, origin='lower', 
                          extent=[mass_unique.min(), mass_unique.max(), ecc_unique.min(), ecc_unique.max()],
                          aspect='auto', interpolation='nearest')
            
            # Add colorbar
            cbar = fig.colorbar(im, ax=ax)
            cbar.set_label(z_label)
            
            # Set labels and title
            ax.set_xlabel('木星质量 (太阳质量倍数)')
            ax.set_ylabel('木星偏心率')
            ax.set_title(title)
            
            # Add grid lines
            ax.set_xticks(mass_unique)
            ax.set_yticks(ecc_unique)
            ax.grid(visible=True, color='w', linestyle='-', linewidth=0.5)
            ax.tick_params(axis='x', rotation=45)
            
            # Update canvas
            canvas.draw()
            
            # Update statistics
            stats_text.delete(1.0, tk.END)
            stats_text.insert(tk.END, f"总共有 {len(simulations)} 次已完成的模拟\n")
            stats_text.insert(tk.END, f"木星质量范围: {mass_unique.min():.2f} 到 {mass_unique.max():.2f}\n")
            stats_text.insert(tk.END, f"木星偏心率范围: {ecc_unique.min():.2f} 到 {ecc_unique.max():.2f}\n")
            stats_text.insert(tk.END, "\n")
            
            if z_axis_var.get() == "柯克伍德带隙清晰度":
                stats_text.insert(tk.END, f"平均带隙清晰度: {gap_clarities.mean():.4f}\n")
                stats_text.insert(tk.END, f"最低带隙清晰度 (最清晰): {gap_clarities.min():.4f}\n")
                stats_text.insert(tk.END, f"最高带隙清晰度 (最模糊): {gap_clarities.max():.4f}\n")
            else:
                stats_text.insert(tk.END, f"平均存活率: {survival_rates.mean():.2f}%\n")
                stats_text.insert(tk.END, f"最高存活率: {survival_rates.max():.2f}%\n")
                stats_text.insert(tk.END, f"最低存活率: {survival_rates.min():.2f}%\n")
            
            conn.close()
            
        except Exception as e:
            messagebox.showerror("错误", f"更新热力图失败: {str(e)}")
    
    refresh_button.config(command=update_heatmap)
    
    # Initialize with first draw
    update_heatmap()
    
    return frame

def main():
    """Main function to run the application"""
    # Create database
    create_database()
    
    # Create main window
    root = tk.Tk()
    root.title("柯克伍德带隙探索器")
    root.geometry("1200x800")
    
    # Create notebook (tabbed interface)
    notebook = ttk.Notebook(root)
    notebook.pack(fill='both', expand=True, padx=5, pady=5)
    
    # Create tabs
    batch_tab = batch_simulator_tab(notebook)
    data_tab = data_mining_tab(notebook)
    
    notebook.add(batch_tab, text="批处理模拟器")
    notebook.add(data_tab, text="数据挖掘与可视化")
    
    # Run main loop
    root.mainloop()

if __name__ == "__main__":
    main()