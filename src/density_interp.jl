using Interpolations
using NPZ

"""
    atmospheric_density(altitude_m, solar_activity=:mean)

Calcola la densità atmosferica in kg/m³ in base all'altitudine e all'attività solare.

# Argomenti
- `altitude_m`: Altitudine in metri (da -2000 a 900000 m)
- `solar_activity`: Livello di attività solare, può essere:
  - `:low` - Bassa attività solare
  - `:mean` - Attività solare media (default)
  - `:high` - Attività solare estremamente alta
  - `:standard` - USA lo standard USSA 1976 (fino a 84852 m)

# Ritorna
- Densità atmosferica in kg/m³

# Esempi
```julia
# Densità al livello del mare
ρ = atmospheric_density(0.0)  # ~1.225 kg/m³

# Densità a 400 km (ISS) con alta attività solare
ρ = atmospheric_density(400000.0, :high)

# Densità all'Everest con standard USSA 1976
ρ = atmospheric_density(8848.0, :standard)
```
"""
function atmospheric_density(altitude_m::Real, solar_activity::Symbol=:mean)
    
    # Verifica che l'attività solare sia valida
    if !(solar_activity in [:low, :mean, :high, :standard])
        error("solar_activity deve essere :low, :mean, :high, o :standard")
    end
    
    # Carica i dati (questo dovrebbe essere fatto una sola volta in pratica)
    data = npzread("atmospheric_data.npz")
    
    # Seleziona i dati appropriati in base all'attività solare
    if solar_activity == :standard
        # Usa U.S. Standard Atmosphere 1976
        altitudes = data["ussa1976_altitude"]
        densities = data["ussa1976_density"]
        max_alt = 84852.0
        min_alt = -2000.0
    elseif solar_activity == :low
        # Usa MSISE-90 bassa attività
        altitudes = data["msise90_low_altitude"]
        densities = data["msise90_low_density"]
        max_alt = 900000.0
        min_alt = 0.0
    elseif solar_activity == :mean
        # Usa MSISE-90 attività media
        altitudes = data["msise90_mean_altitude"]
        densities = data["msise90_mean_density"]
        max_alt = 900000.0
        min_alt = 0.0
    else  # :high
        # Usa MSISE-90 attività alta
        altitudes = data["msise90_high_altitude"]
        densities = data["msise90_high_density"]
        max_alt = 900000.0
        min_alt = 0.0
    end
    
    # Verifica che l'altitudine sia nel range valido
    if altitude_m < min_alt || altitude_m > max_alt
        error("Altitudine fuori dal range valido ($min_alt m - $max_alt m) per solar_activity=$solar_activity")
    end
    
    # Crea interpolatore (lineare o cubico)
    # Per dati atmosferici, l'interpolazione cubica è più accurata
    itp = interpolate((altitudes,), densities, Gridded(Linear()))
    
    # Calcola e ritorna la densità interpolata
    return itp(altitude_m)
end


"""
    atmospheric_density_log(altitude_m, solar_activity=:mean)

Come atmospheric_density, ma usa interpolazione logaritmica.
Più accurata per altitudini elevate dove la densità varia esponenzialmente.

# Esempi
```julia
# Densità a 500 km con interpolazione logaritmica
ρ = atmospheric_density_log(500000.0, :mean)
```
"""
function atmospheric_density_log(altitude_m::Real, solar_activity::Symbol=:mean)
    
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
    
    if altitude_m < min_alt || altitude_m > max_alt
        error("Altitudine fuori dal range valido ($min_alt m - $max_alt m)")
    end
    
    # Interpola il logaritmo della densità (più accurato per variazioni esponenziali)
    log_densities = log.(densities)
    itp = interpolate((altitudes,), log_densities, Gridded(Linear()))
    
    # Ritorna l'esponenziale del valore interpolato
    return exp(itp(altitude_m))
end


"""
    get_atmospheric_properties(altitude_m, solar_activity=:mean)

Ritorna tutte le proprietà atmosferiche disponibili per una data altitudine.

# Ritorna
Un dizionario con:
- `:density` - Densità (kg/m³)
- `:temperature` - Temperatura (K) se disponibile
- `:pressure` - Pressione (Pa) se disponibile
- `:altitude` - Altitudine (m)
- `:solar_activity` - Livello di attività solare
"""
function get_atmospheric_properties(altitude_m::Real, solar_activity::Symbol=:mean)
    
    if !(solar_activity in [:low, :mean, :high, :standard])
        error("solar_activity deve essere :low, :mean, :high, o :standard")
    end
    
    data = npzread("atmospheric_data.npz")
    
    # Seleziona i dati
    if solar_activity == :standard
        altitudes = data["ussa1976_altitude"]
        densities = data["ussa1976_density"]
        temperatures = data["ussa1976_temperature"]
        pressures = data["ussa1976_pressure"]
    elseif solar_activity == :low
        altitudes = data["msise90_low_altitude"]
        densities = data["msise90_low_density"]
        temperatures = data["msise90_low_temperature"]
        pressures = data["msise90_low_pressure"]
    elseif solar_activity == :mean
        altitudes = data["msise90_mean_altitude"]
        densities = data["msise90_mean_density"]
        temperatures = data["msise90_mean_temperature"]
        pressures = data["msise90_mean_pressure"]
    else
        altitudes = data["msise90_high_altitude"]
        densities = data["msise90_high_density"]
        temperatures = data["msise90_high_temperature"]
        pressures = data["msise90_high_pressure"]
    end
    
    # Crea interpolatori
    itp_dens = interpolate((altitudes,), densities, Gridded(Linear()))
    itp_temp = interpolate((altitudes,), temperatures, Gridded(Linear()))
    itp_pres = interpolate((altitudes,), pressures, Gridded(Linear()))
    
    # Ritorna dizionario con tutte le proprietà
    return Dict(
        :altitude => altitude_m,
        :density => itp_dens(altitude_m),
        :temperature => itp_temp(altitude_m),
        :pressure => itp_pres(altitude_m),
        :solar_activity => solar_activity
    )
end


# ============================================================================
# ESEMPI D'USO
# ============================================================================

println("="^70)
println("ESEMPI DI UTILIZZO DELLA FUNZIONE atmospheric_density")
println("="^70)

# Esempio 1: Densità al livello del mare
println("\n1. Densità al livello del mare:")
ρ_sea = atmospheric_density(0.0, :standard)
println("   ρ = $(ρ_sea) kg/m³")

# Esempio 2: Densità all'Everest
println("\n2. Densità al Monte Everest (8,848 m):")
ρ_everest = atmospheric_density(8848.0, :standard)
println("   ρ = $(ρ_everest) kg/m³")

# Esempio 3: Densità alla ISS con diversi livelli di attività solare
println("\n3. Densità alla ISS (~400 km) con diverse attività solari:")
ρ_iss_low = atmospheric_density(400_000.0, :low)
ρ_iss_mean = atmospheric_density(400_000.0, :mean)
ρ_iss_high = atmospheric_density(400_000.0, :high)
println("   Bassa:  ρ = $(ρ_iss_low) kg/m³")
println("   Media:  ρ = $(ρ_iss_mean) kg/m³")
println("   Alta:   ρ = $(ρ_iss_high) kg/m³")
println("   Rapporto high/low: $(ρ_iss_high/ρ_iss_low)")

# Esempio 4: Confronto interpolazione lineare vs logaritmica a 500 km
println("\n4. Confronto interpolazioni a 500 km (attività media):")
ρ_lin = atmospheric_density(500_000.0, :mean)
ρ_log = atmospheric_density_log(500_000.0, :mean)
println("   Lineare:      ρ = $(ρ_lin) kg/m³")
println("   Logaritmica:  ρ = $(ρ_log) kg/m³")
println("   Differenza:   $(abs(ρ_lin - ρ_log)/ρ_log * 100)%")

# Esempio 5: Proprietà complete a 200 km
println("\n5. Proprietà atmosferiche complete a 200 km (attività alta):")
props = get_atmospheric_properties(200_000.0, :high)
println("   Altitudine:   $(props[:altitude]) m")
println("   Densità:      $(props[:density]) kg/m³")
println("   Temperatura:  $(props[:temperature]) K")
println("   Pressione:    $(props[:pressure]) Pa")

# Esempio 6: Profilo di densità
println("\n6. Profilo di densità a varie altitudini (attività media):")
println("   Alt (km)    Densità (kg/m³)")
println("   " * "-"^40)
for alt_km in [0, 10, 50, 100, 200, 400, 600, 800]
    alt_m = alt_km * 1000.0
    if alt_km == 0
        ρ = atmospheric_density(alt_m, :standard)
    else
        ρ = atmospheric_density(alt_m, :mean)
    end
    println("   $(lpad(alt_km, 7))     $(ρ)")
end

println("\n" * "="^70)
println("NOTA: Assicurati che il file 'atmospheric_data.npz' sia")
println("      nella stessa directory del tuo script Julia!")
println("="^70)