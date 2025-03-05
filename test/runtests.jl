using SequentialMeasurements
using Test
using QuantumOptics
using LinearAlgebra
using .SequentialMeasurements.SpectralRepresentation: lf_gaussian, lorentzian

@testset "SequentialMeasurements.jl" begin
    @testset "Core Types" begin
        # Test creation of basic types
        @test begin
            noise_params = SequentialMeasurements.NoiseParams(
                Nb = 0,
                psd_function = lf_gaussian,
                psd_kwargs = Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
                ωu = 2π*1.0,
                Nh = Int(2^10)
            )
            true
        end

        @test begin
            ctrl_params = SequentialMeasurements.CDDParams(
                n = 0,
                t0 = 0.0,
                τg = 1.0,
                nsteps = 10,
                M_max = 5
            )
            true
        end

        @test begin
            measurement_scheme = SequentialMeasurements.MeasurementScheme(
                fout = nothing,
                conditional_evolution = false
            )
            true
        end
    end

    @testset "Spectral Representation" begin
        # Test PSD functions
        @test lorentzian(1.0; σ=1.0, θ=1.0, ω₀=0.0) isa Float64
        @test lorentzian(1.0; σ=2.0, θ=0.5) > 0
        
        @test lf_gaussian(1.0; s=1, Δ0=100, Γ0=0.1) isa Float64
        @test lf_gaussian(1.0; s=1, Δ0=100, Γ0=0.1) > 0
        
        # Test spectral representation with minimal parameters
        t_list = range(0, stop=10, length=100)
        omega_m = range(0.1, stop=10, length=50)
        num_realizations = 2
        
        @test begin
            noise = SequentialMeasurements.SpectralRepresentation.spectral_representation(
                t_list, 
                omega_m, 
                num_realizations, 
                lorentzian;
                σ=1.0, 
                θ=1.0
            )
            size(noise) == (length(t_list), num_realizations)
        end
    end

    @testset "Initialization" begin
        # Test initialization functions
        ctrl_params_dict = Dict(
            :n => 0,
            :t0 => 0.0,
            :τg => 1.0,
            :nsteps => 10,
            :M_max => 5
        )
        
        noise_params_dict = Dict(
            :Nb => 0,
            :psd_function => lf_gaussian,
            :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
            :ωu => 2π*1.0,
            :Nh => Int(2^10)
        )
        
        measurement_scheme_dict = Dict(
            :fout => nothing,
            :conditional_evolution => false
        )
        
        @test begin
            control = SequentialMeasurements.initialize_ctrl(:cdd, ctrl_params_dict)
            control isa SequentialMeasurements.CDDParams
        end
        
        @test begin
            noise = SequentialMeasurements.initialize_noise(:lf_gaussian, noise_params_dict)
            noise isa SequentialMeasurements.NoiseParams
        end
        
        @test begin
            ms = SequentialMeasurements.initialize_measurement_scheme(measurement_scheme_dict)
            ms isa SequentialMeasurements.MeasurementScheme
        end
    end
    
    @testset "CDD Sequence" begin
        # Test the CDD sequence generation
        cdd0 = SequentialMeasurements.ControlledEvolution.cdd(0)
        cdd1 = SequentialMeasurements.ControlledEvolution.cdd(1)
        cdd2 = SequentialMeasurements.ControlledEvolution.cdd(2)
        
        @test cdd0 == [0]
        @test cdd1 == [0, 1, 0, 1]
        @test length(cdd2) > length(cdd1)
        
        # Test that the number of free evolution regions (zeros) in cdd(n) is 2^n
        @test count(x -> x == 0, cdd0) == 2^0
        @test count(x -> x == 0, cdd1) == 2^1
        @test count(x -> x == 0, cdd2) == 2^2
    end
end
