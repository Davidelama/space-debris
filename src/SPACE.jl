using Plots, JSON, Interpolations, NPZ

# Costanti fisiche
const G = 6.67430e-11  # Costante gravitazionale (m³/kg/s²)
const M_terra = 5.972e24  # Massa della Terra (kg)
const R_terra = 6.371e6  # Raggio della Terra (m)

# Funzione per creare file di configurazione di default
function crea_config_default(filename="orbit_config.json")
    config = Dict(
        "simulazione" => Dict(
            "attrito_attivo" => true,
            "tempo_simulazione_giorni" => 30,
            "dt" => 10.0,
            "backend_grafico" => "gr"  # "gr" o "plotly"
        ),
        "orbita" => Dict(
            "altitudine_iniziale_km" => 250.0,
            "inclinazione_gradi" => 51.6,
            "posizione_iniziale" => "automatica"  # o [x, y, z] in km
        ),
        "satellite" => Dict(
            "massa_kg" => 10.0,
            "area_sezione_m2" => 0.1,
            "coefficiente_drag" => 2.2,
            "nome" => "CubeSat-1"
        ),
        "output" => Dict(
            "salva_grafici" => false,
            "prefisso_file" => "orbita",
            "mostra_statistiche" => true,
            "max_punti_grafico" => 5000,  # Massimo numero di punti da plottare
            "salva_dati" => true,  # Salva dati numerici
            "cartella_dati" => "../data",  # Cartella per i file .dat
            "max_punti_dati" => 10000,  # Massimo numero di punti da salvare (0 = tutti)
            "cartella_figure" => "../fig"  # Cartella per i grafici
        )
    )
    
    open(filename, "w") do f
        JSON.print(f, config, 4)
    end
    
    println("✓ File di configurazione creato: $filename")
    return config
end

# Funzione per caricare configurazione
function carica_config(filename="orbit_config.json")
    if !isfile(filename)
        println("⚠️  File $filename non trovato. Creo configurazione di default...")
        return crea_config_default(filename)
    end
    
    config = JSON.parsefile(filename)
    println("✓ Configurazione caricata da: $filename")
    return config
end

# Funzione per calcolare la densità atmosferica in base all'altitudine
function densita_atmosfera_log( solar_activity::Symbol=:mean)
    
    if !(solar_activity in [:low, :mean, :high, :standard])
        error("solar_activity deve essere :low, :mean, :high, o :standard")
    end
    
    data = npzread("atmospheric_data.npz")
    
    if solar_activity == :standard
        altitudes = data["ussa1976_altitude"]
        densities = data["ussa1976_density"]
        max_alt = 84852.0
        min_alt = -2000.0
    elseif solar_activity == :low
        altitudes = data["msise90_low_altitude"]
        densities = data["msise90_low_density"]
        max_alt = 900000.0
        min_alt = 0.0
    elseif solar_activity == :mean
        altitudes = data["msise90_mean_altitude"]
        densities = data["msise90_mean_density"]
        max_alt = 900000.0
        min_alt = 0.0
    else
        altitudes = data["msise90_high_altitude"]
        densities = data["msise90_high_density"]
        max_alt = 900000.0
        min_alt = 0.0
    end

    
    # Interpola il logaritmo della densità (più accurato per variazioni esponenziali)
    log_densities = log.(densities)
    itp = interpolate((altitudes,), log_densities, Gridded(Linear()))
    
    # Ritorna l'esponenziale del valore interpolato
    return itp
end


# Accelerazione totale (gravità + drag)
function accelerazione_totale(pos, vel, params::Dict)
    # Estrai parametri forze dal dizionario
    Cd = params["Cd"]
    A = params["A"]
    m = params["m"]
    attrito_attivo = params["attrito_attivo"]

    # Gravitazione
    r = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)
    a_mag_grav = -G * M_terra / r^2
    a_grav = [a_mag_grav * pos[1] / r,
              a_mag_grav * pos[2] / r,
              a_mag_grav * pos[3] / r]
    
    # Drag atmosferico (se attivo)
    if attrito_attivo
         altitudine = r - R_terra
         rho = exp(atm_intp(altitudine))
         v_mag = sqrt(vel[1]^2 + vel[2]^2 + vel[3]^2)

         if v_mag > 0
             drag_mag = -0.5 * Cd * A * rho * v_mag / m
             a_drag = [drag_mag * vel[1],
                       drag_mag * vel[2],
                       drag_mag * vel[3]]
         else
             a_drag = [0.0, 0.0, 0.0]
         end
     else
         a_drag = [0.0, 0.0, 0.0]
     end

    return a_grav + a_drag
end

# Runge-Kutta 4° ordine
function rk4_step(pos, vel, dt, params::Dict)
    k1_vel = accelerazione_totale(pos, vel, params)
    k1_pos = vel
    
    pos2 = pos + 0.5 * dt * k1_pos
    vel2 = vel + 0.5 * dt * k1_vel
    k2_vel = accelerazione_totale(pos2, vel2, params)
    k2_pos = vel2
    
    pos3 = pos + 0.5 * dt * k2_pos
    vel3 = vel + 0.5 * dt * k2_vel
    k3_vel = accelerazione_totale(pos3, vel3, params)
    k3_pos = vel3
    
    pos4 = pos + dt * k3_pos
    vel4 = vel + dt * k3_vel
    k4_vel = accelerazione_totale(pos4, vel4, params)
    k4_pos = vel4
    
    new_pos = pos + (dt / 6.0) * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos)
    new_vel = vel + (dt / 6.0) * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel)
    
    return new_pos, new_vel
end

# Simulazione principale
function simula_orbita(pos_iniziale, vel_iniziale, tempo_totale, dt, params::Dict)
    n_steps = Int(floor(tempo_totale / dt))
    
    x = zeros(n_steps)
    y = zeros(n_steps)
    z = zeros(n_steps)
    altitudini = zeros(n_steps)
    
    pos = pos_iniziale
    vel = vel_iniziale
    
    for i in 1:n_steps
        x[i] = pos[1]
        y[i] = pos[2]
        z[i] = pos[3]
        
        r = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)
        altitudini[i] = r - R_terra
        
        if altitudini[i] < 100e3
            println("\n⚠️  Satellite rientrato dopo $(round(i*dt/3600, digits=2)) ore!")
            return x[1:i], y[1:i], z[1:i], altitudini[1:i], collect(1:i) * dt
        end
        
        pos, vel = rk4_step(pos, vel, dt, params)
    end
    
    return x, y, z, altitudini, collect(1:n_steps) * dt
end

# ============== MAIN ==============

println("=== Simulatore Orbita Terrestre 3D ===\n")

# Carica configurazione
config = carica_config("orbit_config.json")

# Estrai parametri
sim = config["simulazione"]
orb = config["orbita"]
sat = config["satellite"]
out = config["output"]

# Carica densità atmosferica
atm_intp = densita_atmosfera_log(:mean)

# Imposta backend grafico
if sim["backend_grafico"] == "plotly"
    try
        plotly()
        println("✓ Backend Plotly attivato (interattivo)")
        println("   Usa il mouse per ruotare, zoom e pan nel grafico 3D")
    catch e
        println("⚠️  Plotly non disponibile: $e")
        println("   Installa con: using Pkg; Pkg.add(\"PlotlyJS\")")
        println("   Uso backend GR invece...")
        gr()
    end
else
    gr()
    println("✓ Backend GR attivato (statico)")
end

# Calcola parametri orbitali
altitudine_iniziale = orb["altitudine_iniziale_km"] * 1e3
r_orbita = R_terra + altitudine_iniziale
inclinazione = orb["inclinazione_gradi"] * π / 180

pos_iniziale = [r_orbita, 0.0, 0.0]
v_orbitale = sqrt(G * M_terra / r_orbita)
vel_iniziale = [0.0, 
                v_orbitale * cos(inclinazione), 
                v_orbitale * sin(inclinazione)]

# Parametri satellite
massa = sat["massa_kg"]
area = sat["area_sezione_m2"]
Cd = sat["coefficiente_drag"]
attrito = sim["attrito_attivo"]

# Dizionario parametri forze (passato alle funzioni di integrazione) — ora include Cd
params_forze = Dict("A" => area, "m" => massa, "attrito_attivo" => attrito, "Cd" => Cd)

# Tempo simulazione
periodo = 2 * π * r_orbita / v_orbitale
tempo_totale = sim["tempo_simulazione_giorni"] * 24 * 3600
dt = sim["dt"]

# Mostra info
if out["mostra_statistiche"]
    println("\n--- Parametri Simulazione ---")
    println("Satellite: $(sat["nome"])")
    println("Altitudine iniziale: $(orb["altitudine_iniziale_km"]) km")
    println("Inclinazione: $(orb["inclinazione_gradi"])°")
    println("Velocità orbitale: $(round(v_orbitale/1000, digits=2)) km/s")
    println("Periodo: $(round(periodo/60, digits=1)) min")
    println("Massa: $(massa) kg")
    println("Area: $(area) m²")
    println("Drag: $(attrito ? "ATTIVO" : "DISATTIVO")")
    println("Tempo simulazione: $(sim["tempo_simulazione_giorni"]) giorni")
    println("\nSimulazione in corso...")
end

# Esegui simulazione (passare il dizionario dei parametri forze, che contiene Cd)
x, y, z, altitudini, t = simula_orbita(pos_iniziale, vel_iniziale, tempo_totale, dt,
                                     params_forze)

# Salva dati se richiesto
if out["salva_dati"]
    println("\nSalvataggio dati...")
    
    # Downsampling per salvataggio dati
    n_punti_totali = length(x)
    max_punti_dati = out["max_punti_dati"]
    
    if max_punti_dati > 0 && n_punti_totali > max_punti_dati
        step_dati = Int(ceil(n_punti_totali / max_punti_dati))
        x_save = x[1:step_dati:end]
        y_save = y[1:step_dati:end]
        z_save = z[1:step_dati:end]
        alt_save = altitudini[1:step_dati:end]
        dt_save = dt * step_dati
        println("  Downsampling dati: usando 1 punto ogni $step_dati ($(length(x_save)) punti)")
    else
        x_save = x
        y_save = y
        z_save = z
        alt_save = altitudini
        dt_save = dt
        println("  Salvando tutti i $(n_punti_totali) punti")
    end
    
    # Crea nome cartella con parametri
    params_str = "alt$(Int(round(orb["altitudine_iniziale_km"])))km_" *
                 "inc$(Int(round(orb["inclinazione_gradi"])))deg_" *
                 "m$(massa)kg_" *
                 "A$(area)m2_" *
                 "Cd$(Cd)_" *
                 "$(attrito ? "drag" : "nodrag")"
    
    cartella_output = joinpath(out["cartella_dati"], params_str)
    
    # Crea cartella se non esiste
    if !isdir(cartella_output)
        mkpath(cartella_output)
        println("  ✓ Cartella creata: $cartella_output")
    end
    
    # Salva posizioni (x, y, z in km)
    open(joinpath(cartella_output, "posizione_$params_str.dat"), "w") do f
        println(f, "# Posizione del satellite")
        println(f, "# Colonne: tempo(s)  x(km)  y(km)  z(km)")
        for i in 1:length(x_save)
            println(f, "$(i*dt_save)  $(x_save[i]/1e3)  $(y_save[i]/1e3)  $(z_save[i]/1e3)")
        end
    end
    
    # Salva altitudine
    open(joinpath(cartella_output, "altitudine_$params_str.dat"), "w") do f
        println(f, "# Altitudine del satellite")
        println(f, "# Colonne: tempo(s)  altitudine(km)")
        for i in 1:length(alt_save)
            println(f, "$(i*dt_save)  $(alt_save[i]/1e3)")
        end
    end
    
    # Salva parametri di configurazione
    open(joinpath(cartella_output, "parametri_$params_str.txt"), "w") do f
        println(f, "=== Parametri Simulazione ===")
        println(f, "Satellite: $(sat["nome"])")
        println(f, "Altitudine iniziale: $(orb["altitudine_iniziale_km"]) km")
        println(f, "Inclinazione: $(orb["inclinazione_gradi"])°")
        println(f, "Massa: $(massa) kg")
        println(f, "Area sezione: $(area) m²")
        println(f, "Coefficiente drag: $(Cd)")
        println(f, "Drag atmosferico: $(attrito ? "ATTIVO" : "DISATTIVO")")
        println(f, "Tempo simulazione: $(sim["tempo_simulazione_giorni"]) giorni")
        println(f, "Passo temporale calcolo: $(dt) s")
        println(f, "Passo temporale salvato: $(dt_save) s")
        println(f, "Numero punti calcolati: $(n_punti_totali)")
        println(f, "Numero punti salvati: $(length(x_save))")
        if attrito
            println(f, "\n=== Risultati Decadimento ===")
            println(f, "Altitudine finale: $(round(altitudini[end]/1e3, digits=1)) km")
            println(f, "Perdita totale: $(round((altitudini[1]-altitudini[end])/1e3, digits=1)) km")
        end
    end
    
    println("  ✓ Dati salvati in: $cartella_output")
    println("    - posizione_$params_str.dat")
    println("    - altitudine_$params_str.dat")
    println("    - parametri_$params_str.txt")
end

# Visualizzazione 3D
println("Creazione grafici...")

# Downsampling per migliorare performance grafiche
n_punti = length(x)
max_punti = out["max_punti_grafico"]
step = max(1, Int(ceil(n_punti / max_punti)))

if step > 1
    println("  Downsampling: usando 1 punto ogni $step ($(Int(floor(n_punti/step))) punti totali)")
end

x_plot = x[1:step:end]
y_plot = y[1:step:end]
z_plot = z[1:step:end]
alt_plot = altitudini[1:step:end]
tempo_plot = collect(1:length(alt_plot)) * dt * step / 3600

# Terra
u = range(0, 2π, length=50)
v = range(0, π, length=50)
x_terra = R_terra * cos.(u) .* sin.(v)'
y_terra = R_terra * sin.(u) .* sin.(v)'
z_terra = R_terra * ones(length(u)) .* cos.(v)'

plt1 = plot(x_plot/1e6, y_plot/1e6, z_plot/1e6,
    label="Traiettoria $(sat["nome"])",
    linewidth=2,
    color=:blue,
    xlabel="x (1000 km)",
    ylabel="y (1000 km)",
    zlabel="z (1000 km)",
    title="Orbita 3D - $(attrito ? "Con" : "Senza") Drag",
    legend=:topright,
    size=(900, 800),
    camera=(30, 30))

surface!(x_terra/1e6, y_terra/1e6, z_terra/1e6,
    alpha=0.3,
    color=:lightblue,
    label="Terra",
    colorbar=false)

scatter!([x_plot[1]/1e6], [y_plot[1]/1e6], [z_plot[1]/1e6],
    label="Inizio",
    markersize=6,
    color=:green)

if !out["salva_grafici"]
    display(plt1)
end

# Grafico altitudine
if attrito
    plt2 = plot(tempo_plot, alt_plot/1e3,
        xlabel="Tempo (ore)",
        ylabel="Altitudine (km)",
        title="Decadimento Orbitale - $(sat["nome"])",
        linewidth=2,
        color=:red,
        legend=false,
        size=(900, 400))
    
    hline!([100], linestyle=:dash, color=:black, linewidth=2)
    
    if !out["salva_grafici"]
        display(plt2)
    end
    
    if out["mostra_statistiche"]
        alt_i = altitudini[1] / 1e3
        alt_f = altitudini[end] / 1e3
        perdita = alt_i - alt_f
        
        println("\n=== Statistiche Decadimento ===")
        println("Altitudine iniziale: $(round(alt_i, digits=1)) km")
        println("Altitudine finale: $(round(alt_f, digits=1)) km")
        println("Perdita totale: $(round(perdita, digits=1)) km")
        println("Perdita media: $(round(perdita/last(t)*3600*24, digits=2)) km/giorno")
    end
    
end

if out["salva_grafici"]
    println("\nSalvataggio grafici...")
    
    # Crea nome cartella con parametri (stesso della cartella dati)
    params_str = "alt$(Int(round(orb["altitudine_iniziale_km"])))km_" *
                 "inc$(Int(round(orb["inclinazione_gradi"])))deg_" *
                 "m$(massa)kg_" *
                 "A$(area)m2_" *
                 "Cd$(Cd)_" *
                 "$(attrito ? "drag" : "nodrag")"
    
    cartella_figure = joinpath(out["cartella_figure"], params_str)
    
    # Crea cartella se non esiste
    if !isdir(cartella_figure)
        mkpath(cartella_figure)
        println("  ✓ Cartella figure creata: $cartella_figure")
    end
    
    # Salva grafico 3D
    savefig(plt1, joinpath(cartella_figure, "traiettoria3d_$params_str.png"))
    println("  ✓ Grafico 3D salvato: traiettoria3d_$params_str.png")
    
    # Salva grafico decadimento se drag attivo
    if attrito
        savefig(plt2, joinpath(cartella_figure, "decadimento_$params_str.png"))
        println("  ✓ Grafico decadimento salvato: decadimento_$params_str.png")
    end
else
    println("\n(I grafici sono mostrati a schermo, non salvati)")
end

println("\n✓ Simulazione completata!")
println("Premi INVIO per chiudere...")
readline()