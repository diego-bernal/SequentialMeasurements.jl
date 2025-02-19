using Test
using QuantumOptics
using SequentialMeasurements
include("../scripts/error_virtualization.jl")

@testset "Error Virtualization Tests" begin
    # Test initialization of parameters
    @testset "Parameter Initialization" begin
        ctrl_params_dict = Dict(
            :n => 0,
            :t0 => 0.0,
            :τg => 1.0,
            :nsteps => 10,
            :M_max => 60
        )
        
        noise_params_dict = Dict(
            :Nb => 0,
            :psd_function => lf_gaussian,
            :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
            :ωu => 2π*1.0,
            :Nh => Int(2.0^12)
        )
        
        measurement_scheme_dict = Dict(
            :fout => process_fidelity_with_plus,
            :conditional_evolution => false
        )
        
        @test_nowarn control = initialize_ctrl(:cdd, ctrl_params_dict)
        @test_nowarn noise = initialize_noise(:lf_gaussian, noise_params_dict)
        @test_nowarn measurement_scheme = initialize_measurement_scheme(measurement_scheme_dict)
    end

    # Test spectral density functions
    @testset "Spectral Functions" begin
        ω = range(0.1, stop=10.0, length=100)
        @test all(lf_gaussian.(ω; s=1.0, Δ0=50, Γ0=0.1) .>= 0)
    end

    # Test noise generation
    @testset "Noise Generation" begin
        n_runs = 2
        Nb = 0
        ωu = 2π*1.0
        Nh = Int(2.0^8)  # Reduced for testing
        
        result = noise_generation_loop(
            n_runs, 
            Nb,
            lf_gaussian,
            ωu, 
            Nh;
            s=1.0, 
            Δ0=50, 
            Γ0=0.1
        )
        
        @test length(result[1]) == n_runs  # Test ξ_sys length
        @test isempty(result[2])  # Test ξ_bath is empty
    end

    # Test single realization
    @testset "Single Realization" begin
        # Initialize parameters
        control = initialize_ctrl(:cdd, Dict(
            :n => 0,
            :t0 => 0.0,
            :τg => 1.0,
            :nsteps => 10,
            :M_max => 2
        ))
        
        noise = initialize_noise(:lf_gaussian, Dict(
            :Nb => 0,
            :psd_function => lf_gaussian,
            :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
            :ωu => 2π*1.0,
            :Nh => Int(2.0^8)
        ))
        
        measurement_scheme = initialize_measurement_scheme(Dict(
            :fout => process_fidelity_with_plus,
            :conditional_evolution => false
        ))
        
        # Run a single realization
        num_realizations = 1
        result = run_simulation(
            num_realizations,
            noise,
            control,
            tempdir(),  # Use temporary directory for testing
            measurement_scheme,
            error_virtualization_realization;
            save_output=false
        )
        
        @test length(result) > 0  # Should return some results
    end
end
