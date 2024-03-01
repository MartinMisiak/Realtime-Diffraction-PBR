#ifndef UNIVERSAL_BRDF_INCLUDED
#define UNIVERSAL_BRDF_INCLUDED

#include "Packages/com.unity.render-pipelines.core/ShaderLibrary/BSDF.hlsl"
#include "Packages/com.unity.render-pipelines.core/ShaderLibrary/CommonMaterial.hlsl"
#include "Packages/com.unity.render-pipelines.core/ShaderLibrary/Random.hlsl"
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Deprecated.hlsl"
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/SurfaceData.hlsl"


#ifdef _DIFFRACTION_PATTERN_OPEN_SIMPLEX_2
    #include "Assets/Shaders/OpenSimplex2.hlsl"
#else
    #include "Assets/Shaders/SimplexNoise.hlsl"
#endif




#define kDielectricSpec half4(0.04, 0.04, 0.04, 1.0 - 0.04) // standard dielectric reflectivity coef at incident angle (= 4%)

struct BRDFData
{
    half3 albedo;
    half3 diffuse;
    half3 specular;
    half reflectivity;
    half perceptualRoughness;
    half roughness;
    half roughness2;
    half grazingTerm;

    // We save some light invariant BRDF terms so we don't have to recompute
    // them in the light loop. Take a look at DirectBRDF function for detailed explaination.
    half normalizationTerm;     // roughness * 4.0 + 2.0
    half roughness2MinusOne;    // roughness^2 - 1.0
    
    // Clausen BRDF
#if defined _DIRECT_LIGHT_BRDF_DIFFRACTION
    float2 uvs;
    float3 wsPosition;
#endif
    // Clausen BRDF
    
};

half ReflectivitySpecular(half3 specular)
{
#if defined(SHADER_API_GLES)
    return specular.r; // Red channel - because most metals are either monocrhome or with redish/yellowish tint
#else
    return Max3(specular.r, specular.g, specular.b);
#endif
}

half OneMinusReflectivityMetallic(half metallic)
{
    // We'll need oneMinusReflectivity, so
    //   1-reflectivity = 1-lerp(dielectricSpec, 1, metallic) = lerp(1-dielectricSpec, 0, metallic)
    // store (1-dielectricSpec) in kDielectricSpec.a, then
    //   1-reflectivity = lerp(alpha, 0, metallic) = alpha + metallic*(0 - alpha) =
    //                  = alpha - metallic * alpha
    half oneMinusDielectricSpec = kDielectricSpec.a;
    return oneMinusDielectricSpec - metallic * oneMinusDielectricSpec;
}

half MetallicFromReflectivity(half reflectivity)
{
    half oneMinusDielectricSpec = kDielectricSpec.a;
    return (reflectivity - kDielectricSpec.r) / oneMinusDielectricSpec;
}

inline void InitializeBRDFDataDirect(half3 albedo, half3 diffuse, half3 specular, half reflectivity, half oneMinusReflectivity, half smoothness, inout half alpha, out BRDFData outBRDFData)
{
    outBRDFData = (BRDFData)0;
    outBRDFData.albedo = albedo;
    outBRDFData.diffuse = diffuse;
    outBRDFData.specular = specular;
    outBRDFData.reflectivity = reflectivity;

    outBRDFData.perceptualRoughness = PerceptualSmoothnessToPerceptualRoughness(smoothness);
    
    #if defined _DIRECT_LIGHT_BRDF_DIFFRACTION || _DIRECT_LIGHT_BRDF_COOK_TORRANCE
        outBRDFData.roughness       = outBRDFData.perceptualRoughness;
    #else
        outBRDFData.roughness       = max(PerceptualRoughnessToRoughness(outBRDFData.perceptualRoughness), HALF_MIN_SQRT);
    #endif
    
    outBRDFData.roughness2          = max(outBRDFData.roughness * outBRDFData.roughness, HALF_MIN);
    outBRDFData.grazingTerm         = saturate(smoothness + reflectivity);
    outBRDFData.normalizationTerm   = outBRDFData.roughness * half(4.0) + half(2.0);
    outBRDFData.roughness2MinusOne  = outBRDFData.roughness2 - half(1.0);

#ifdef _ALPHAPREMULTIPLY_ON
    outBRDFData.diffuse *= alpha;
    alpha = alpha * oneMinusReflectivity + reflectivity; // NOTE: alpha modified and propagated up.
#endif
}

// Legacy: do not call, will not correctly initialize albedo property.
inline void InitializeBRDFDataDirect(half3 diffuse, half3 specular, half reflectivity, half oneMinusReflectivity, half smoothness, inout half alpha, out BRDFData outBRDFData)
{
    InitializeBRDFDataDirect(half3(0.0, 0.0, 0.0), diffuse, specular, reflectivity, oneMinusReflectivity, smoothness, alpha, outBRDFData);
}

// Initialize BRDFData for material, managing both specular and metallic setup using shader keyword _SPECULAR_SETUP.
inline void InitializeBRDFData(half3 albedo, half metallic, half3 specular, half smoothness, inout half alpha, out BRDFData outBRDFData)
{
#ifdef _SPECULAR_SETUP
    half reflectivity = ReflectivitySpecular(specular);
    half oneMinusReflectivity = half(1.0) - reflectivity;
    half3 brdfDiffuse = albedo * (half3(1.0, 1.0, 1.0) - specular);
    half3 brdfSpecular = specular;
#else
    half oneMinusReflectivity = OneMinusReflectivityMetallic(metallic);
    half reflectivity = half(1.0) - oneMinusReflectivity;
    half3 brdfDiffuse = albedo * oneMinusReflectivity;
    half3 brdfSpecular = lerp(kDieletricSpec.rgb, albedo, metallic);
#endif

    InitializeBRDFDataDirect(albedo, brdfDiffuse, brdfSpecular, reflectivity, oneMinusReflectivity, smoothness, alpha, outBRDFData);
}

inline void InitializeBRDFData(inout SurfaceData surfaceData, out BRDFData brdfData)
{
    InitializeBRDFData(surfaceData.albedo, surfaceData.metallic, surfaceData.specular, surfaceData.smoothness, surfaceData.alpha, brdfData);
}

half3 ConvertF0ForClearCoat15(half3 f0)
{
#if defined(SHADER_API_MOBILE)
    return ConvertF0ForAirInterfaceToF0ForClearCoat15Fast(f0);
#else
    return ConvertF0ForAirInterfaceToF0ForClearCoat15(f0);
#endif
}

inline void InitializeBRDFDataClearCoat(half clearCoatMask, half clearCoatSmoothness, inout BRDFData baseBRDFData, out BRDFData outBRDFData)
{
    outBRDFData = (BRDFData)0;
    outBRDFData.albedo = half(1.0);

    // Calculate Roughness of Clear Coat layer
    outBRDFData.diffuse             = kDielectricSpec.aaa; // 1 - kDielectricSpec
    outBRDFData.specular            = kDielectricSpec.rgb;
    outBRDFData.reflectivity        = kDielectricSpec.r;

    outBRDFData.perceptualRoughness = PerceptualSmoothnessToPerceptualRoughness(clearCoatSmoothness);
    outBRDFData.roughness           = max(PerceptualRoughnessToRoughness(outBRDFData.perceptualRoughness), HALF_MIN_SQRT);
    outBRDFData.roughness2          = max(outBRDFData.roughness * outBRDFData.roughness, HALF_MIN);
    outBRDFData.normalizationTerm   = outBRDFData.roughness * half(4.0) + half(2.0);
    outBRDFData.roughness2MinusOne  = outBRDFData.roughness2 - half(1.0);
    outBRDFData.grazingTerm         = saturate(clearCoatSmoothness + kDielectricSpec.x);

// Relatively small effect, cut it for lower quality
#if !defined(SHADER_API_MOBILE)
    // Modify Roughness of base layer using coat IOR
    half ieta                        = lerp(1.0h, CLEAR_COAT_IETA, clearCoatMask);
    half coatRoughnessScale          = Sq(ieta);
    half sigma                       = RoughnessToVariance(PerceptualRoughnessToRoughness(baseBRDFData.perceptualRoughness));

    baseBRDFData.perceptualRoughness = RoughnessToPerceptualRoughness(VarianceToRoughness(sigma * coatRoughnessScale));

    // Recompute base material for new roughness, previous computation should be eliminated by the compiler (as it's unused)
    baseBRDFData.roughness          = max(PerceptualRoughnessToRoughness(baseBRDFData.perceptualRoughness), HALF_MIN_SQRT);
    baseBRDFData.roughness2         = max(baseBRDFData.roughness * baseBRDFData.roughness, HALF_MIN);
    baseBRDFData.normalizationTerm  = baseBRDFData.roughness * 4.0h + 2.0h;
    baseBRDFData.roughness2MinusOne = baseBRDFData.roughness2 - 1.0h;
#endif

    // Darken/saturate base layer using coat to surface reflectance (vs. air to surface)
    baseBRDFData.specular = lerp(baseBRDFData.specular, ConvertF0ForClearCoat15(baseBRDFData.specular), clearCoatMask);
    // TODO: what about diffuse? at least in specular workflow diffuse should be recalculated as it directly depends on it.
}

BRDFData CreateClearCoatBRDFData(SurfaceData surfaceData, inout BRDFData brdfData)
{
    BRDFData brdfDataClearCoat = (BRDFData)0;

    #if defined(_CLEARCOAT) || defined(_CLEARCOATMAP)
    // base brdfData is modified here, rely on the compiler to eliminate dead computation by InitializeBRDFData()
    InitializeBRDFDataClearCoat(surfaceData.clearCoatMask, surfaceData.clearCoatSmoothness, brdfData, brdfDataClearCoat);
    #endif

    return brdfDataClearCoat;
}

// Computes the specular term for EnvironmentBRDF
half3 EnvironmentBRDFSpecular(BRDFData brdfData, half fresnelTerm)
{
    float surfaceReduction = 1.0 / (brdfData.roughness2 + 1.0);
    return half3(surfaceReduction * lerp(brdfData.specular, brdfData.grazingTerm, fresnelTerm));
}

half3 EnvironmentBRDF(BRDFData brdfData, half3 indirectDiffuse, half3 indirectSpecular, half fresnelTerm)
{
    half3 c = indirectDiffuse * brdfData.diffuse;
    c += indirectSpecular * EnvironmentBRDFSpecular(brdfData, fresnelTerm);
    return c;
}

// Environment BRDF without diffuse for clear coat
half3 EnvironmentBRDFClearCoat(BRDFData brdfData, half clearCoatMask, half3 indirectSpecular, half fresnelTerm)
{
    float surfaceReduction = 1.0 / (brdfData.roughness2 + 1.0);
    return indirectSpecular * EnvironmentBRDFSpecular(brdfData, fresnelTerm) * clearCoatMask;
}

// Computes the scalar specular term for Minimalist CookTorrance BRDF
// NOTE: needs to be multiplied with reflectance f0, i.e. specular color to complete
half DirectBRDFSpecular(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS)
{
    float3 lightDirectionWSFloat3 = float3(lightDirectionWS);
    float3 halfDir = SafeNormalize(lightDirectionWSFloat3 + float3(viewDirectionWS));

    float NoH = saturate(dot(float3(normalWS), halfDir));
    half LoH = half(saturate(dot(lightDirectionWSFloat3, halfDir)));

    // GGX Distribution multiplied by combined approximation of Visibility and Fresnel
    // BRDFspec = (D * V * F) / 4.0
    // D = roughness^2 / ( NoH^2 * (roughness^2 - 1) + 1 )^2
    // V * F = 1.0 / ( LoH^2 * (roughness + 0.5) )
    // See "Optimizing PBR for Mobile" from Siggraph 2015 moving mobile graphics course
    // https://community.arm.com/events/1155

    // Final BRDFspec = roughness^2 / ( NoH^2 * (roughness^2 - 1) + 1 )^2 * (LoH^2 * (roughness + 0.5) * 4.0)
    // We further optimize a few light invariant terms
    // brdfData.normalizationTerm = (roughness + 0.5) * 4.0 rewritten as roughness * 4.0 + 2.0 to a fit a MAD.
    float d = NoH * NoH * brdfData.roughness2MinusOne + 1.00001f;
    half d2 = half(d * d);

    half LoH2 = LoH * LoH;
    half specularTerm = brdfData.roughness2 / (d2 * max(half(0.1), LoH2) * brdfData.normalizationTerm);

    // On platforms where half actually means something, the denominator has a risk of overflow
    // clamp below was added specifically to "fix" that, but dx compiler (we convert bytecode to metal/gles)
    // sees that specularTerm have only non-negative terms, so it skips max(0,..) in clamp (leaving only min(100,...))
#if defined (SHADER_API_MOBILE) || defined (SHADER_API_SWITCH)
    specularTerm = specularTerm - HALF_MIN;
    specularTerm = clamp(specularTerm, 0.0, 100.0); // Prevent FP16 overflow on mobiles
#endif

return specularTerm;
}


// Based on Minimalist CookTorrance BRDF
// Implementation is slightly different from original derivation: http://www.thetenthplanet.de/archives/255
//
// * NDF [Modified] GGX
// * Modified Kelemen and Szirmay-Kalos for Visibility term
// * Fresnel approximated with 1/LdotH
half3 DirectBDRF(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS, bool specularHighlightsOff)
{
    // Can still do compile-time optimisation.
    // If no compile-time optimized, extra overhead if branch taken is around +2.5% on some untethered platforms, -10% if not taken.
    [branch] if (!specularHighlightsOff)
    {
        half specularTerm = DirectBRDFSpecular(brdfData, normalWS, lightDirectionWS, viewDirectionWS);
        half3 color = brdfData.diffuse + specularTerm * brdfData.specular;
        return color;
    }
    else
        return brdfData.diffuse;
}

// Based on Minimalist CookTorrance BRDF
// Implementation is slightly different from original derivation: http://www.thetenthplanet.de/archives/255
//
// * NDF [Modified] GGX
// * Modified Kelemen and Szirmay-Kalos for Visibility term
// * Fresnel approximated with 1/LdotH
half3 DirectBRDF(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS)
{
#ifndef _SPECULARHIGHLIGHTS_OFF
    return brdfData.diffuse + DirectBRDFSpecular(brdfData, normalWS, lightDirectionWS, viewDirectionWS) * brdfData.specular;
#else
    return brdfData.diffuse;
#endif
}


// Custom BRDF-Code
// ------------------------------------------------------

// Performs a Cholesky-Decomposition of inputMat and returns the lower matrix
// Source: https://www.lume.ufrgs.br/bitstream/handle/10183/151001/001009773.pdf
float3x3 CholeskyBanachiewicz3x3Decomp(float3x3 inputMat)
{
    float3x3 L =
    {
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    };
    for (int i = 0; i < 3; i++)
    {
    
        for (int j = 0; j < (i + 1); j++)
        {
            float sum = 0;
            for (int k = 0; k < j; k++)
                sum += L[i][k] * L[j][k];

            if (i == j)
                L[i][j] = sqrt(inputMat[i][i] - sum);
            else
                L[i][j] = ((1.0 / L[j][j]) * (inputMat[i][j] - sum));
        }
    }
    return L;
}

float3 cartesian2Polar(float3 cartPos)
{
    float radius = length(cartPos);
    float theta  = atan2(cartPos.y, cartPos.x);
    float phi    = acos(cartPos.z / radius);
    return float3(theta, phi, radius);
}

float3 fresnelSchlick(float3 F0, float vFresnel, float VdotH)
{
    return (F0 + (1.0 - F0) * pow(1.0 - clamp(VdotH, 0.0, 1.0), 5.0));
}

// Based on:
// Fresnel term approximations for metals [Lazanzi et al. 2005]
// Fresnel Equations Considered Harmful [Hoffman 2019]
// TODO: Investigate. Formula has a very specific parametrization...not general
float3 fresnelSchlickLazanyi(float3 F0, float vFresnel, float VdotH)
{
    float3 a = float3(1.39119566, 1.0848784, 0.71264977); // Corresponds to alluminium and some Olaf magic
    return fresnelSchlick(F0, vFresnel, VdotH) - a * VdotH * pow(1 - VdotH, 6.0);
}

float distributionGTR2(float NdotH, float a)
{
    float a2 = a * a;
    float t = 1 + (a2 - 1) * NdotH * NdotH;
    return a2 / (PI * t * t);
}

float geometricSmithGGX(float NdotX, float alphaG)
{
    float a = alphaG * alphaG;
    float b = NdotX * NdotX;
    return 1.0 / (NdotX + sqrt(a + b - a * b));
}

// Based on [Clausen et al. 2022]
float3 shift_function(float NdotH, float w, float h)
{
    float theta_m = acos(NdotH);
    float m = h * cos(w * theta_m);
//  srgb values
    float m_r = 42.45318742;
    float m_g = -56.98651893;
    float m_b = -159.23105974;
    // adobe rgb
//    float m_r = 14.12228819;
//    float m_g = -56.99255935;
//    float m_b = -155.01640388;

    float3 shift_rgb = float3(m_r * m + 1.0, m_g * m + 1.0, m_b * m + 1.0);
    return shift_rgb;
}

float cov_model(float NdotH, float w)
{
    float theta_m = acos(NdotH);
    float slope = cos(w * theta_m);
    slope = (slope + 1.0) * 0.5;
    
    return slope;
}

// Specular, Cook-Torrance BRDF used in Clausen et al. without Diffraction pattern
float3 DirectBRDFSpecular_GGXReference(BRDFData brdfData, float3 normalWS, float3 lightDirectionWS, float3 viewDirectionWS)
{
    float3 halfDir = SafeNormalize(lightDirectionWS + viewDirectionWS);
    float VdotH = dot(viewDirectionWS, halfDir); // saturate ?
    float LdotH = saturate(dot(lightDirectionWS, halfDir));
    float NdotV = dot(normalWS, viewDirectionWS); // saturate ?
    float NdotL = dot(normalWS, lightDirectionWS); // saturate ?
    float NdotH = saturate(dot(float3(normalWS), halfDir));
    float3 F0       = brdfData.specular;
    float roughness = brdfData.roughness;
    
    //float3 F = fresnelSchlickLazanyi(F0, 0, VdotH);  // Why should we use Lazanyis improvement to metals ?
    float3 F = fresnelSchlick(F0, 0, VdotH);
    float Ds = distributionGTR2(NdotH, roughness);
    float Gs = geometricSmithGGX(NdotL, roughness) * geometricSmithGGX(NdotV, roughness);
    float3 BRDFspecular = F * Ds * Gs;
    
    return BRDFspecular;
}

// Also performs run-time filtering via multisampling(4x)
// For stereoscopic rendering, halfVecWs should be computed using the cyclopean eye position
float3 sampleDiffractionPattern(float2 uvs, float3 halfVecWs, float3 normalWs)
{
    const float UV_TO_SPECKLE_FACTOR = 0.5; // How large is a speckle in the used noise function (speckle size = 1/2 period)
    const float M = _DiffractionZW_Scale;
    
    // If more than 1 speckle fits into the dimension of a screen pixel, the speckle intensity is reduced linearly to their density.
    // This value describes how many additional speckles can be present in the dimension of a screen pixel, before reaching 0 intensity.
    // We use a value of 1 here, as we can fit a total of 2 speckles per pixel dimension without introducing aliasing (due to 4x multisampling)
    // however due to the low contrast of the speckle pattern, we could go up to a value of 2 without any noticeable aliasing
    const float AMPLITUDE_REDUCTION_FALLOFF = 1.45; // MAX_SPECKLES_PER_PIXEL = (AMPLITUDE_REDUCTION_FACTOR + 1)^2
    
    // Compute screen-space derivative of uvs -> How much do uvs change from one pixel to the next ?
    float2 dx_vtc = ddx(uvs);
    float2 dy_vtc = ddy(uvs);
    float delta_max_sqr = max(dot(dx_vtc, dx_vtc), dot(dy_vtc, dy_vtc));
    float delta_uv = sqrt(delta_max_sqr);
    // How many speckles fit along one dimension of a screen pixel
    float sqrt_speckles_per_pixel = delta_uv / UV_TO_SPECKLE_FACTOR;
    // 0 to 1 speckles_per_pixel -> no modulation. 1 to (1 + AMP_RED_FALLOFF) -> linear reduction
    float amplitude_modulation = 1.0 - min((max(sqrt_speckles_per_pixel - 1, 0) / AMPLITUDE_REDUCTION_FALLOFF), 1.0);
    

    // Determine the last two dimensions of the 4D lookup via polar coordinate deltas
    float3 h_polar = cartesian2Polar(halfVecWs);
    float3 n_polar = cartesian2Polar(normalWs);
    float h_a = h_polar.x - n_polar.x;
    float h_p = h_polar.y - n_polar.y;
    h_a *= M;
    h_p *= M;
#if defined(USING_STEREO_MATRICES)
        // addtionally offset polar angle by a fixed value for each eye. This should result in a comfortable highlight disparity
        h_a += _DiffractionStereoSpecularity * UV_TO_SPECKLE_FACTOR * unity_StereoEyeIndex;
#endif
    
        
    float3 filtered_noise = float3(0, 0, 0);
    if (sqrt_speckles_per_pixel > 1.0 && (sqrt_speckles_per_pixel <= (1 + AMPLITUDE_REDUCTION_FALLOFF))) // 4x multi-sampling
    {
        for (int i = 0; i < 4; i++)
        {
            float2 sample_location = uvs + MASK_2X2_GRID[i] * delta_uv;
            float3 sampled_noise = float3(snoise(float4(sample_location, float2(h_a, h_p))),
                                            snoise(float4(sample_location + 43, float2(h_a, h_p))),
                                            snoise(float4(sample_location - 17, float2(h_a, h_p))));
        
            filtered_noise += sampled_noise;
        }
        filtered_noise /= 4.0;
    }
    else if (sqrt_speckles_per_pixel <= 1) // no multi-sampling
    {
        filtered_noise = float3(snoise(float4(uvs, float2(h_a, h_p))),
                                snoise(float4(uvs + 43, float2(h_a, h_p))),
                                snoise(float4(uvs - 17, float2(h_a, h_p))));
    }
    
    // Modulate noise amplitude towards 0, to converge against mean
    return amplitude_modulation * filtered_noise;
}

#if defined _DIRECT_LIGHT_BRDF_DIFFRACTION
// Specular, Cook-Torrance BRDF used in Clausen et al. with Diffraction pattern
float3 DirectBRDFSpecular_Diffraction(BRDFData brdfData, float3 normalWS, float3 lightDirectionWS, float3 viewDirectionWS)
{
    float3 halfDir = SafeNormalize(lightDirectionWS + viewDirectionWS);
    float VdotH = dot(viewDirectionWS, halfDir); // saturate ?
    float LdotH = saturate(dot(lightDirectionWS, halfDir));
    float NdotV = dot(normalWS, viewDirectionWS); // saturate ?
    float NdotL = dot(normalWS, lightDirectionWS); // saturate ?
    float NdotH = saturate(dot(float3(normalWS), halfDir));
    float3 F0   = brdfData.specular;
    float roughness = brdfData.roughness;
    
    float3 shift      = shift_function(NdotH, _DiffractionWidth, _DiffractionHeight); //shift = float3(0.40367881, 0.3801785, 0.34816286);
    float3 rnd_number = 0;

    #if defined _DIRECT_LIGHT_BRDF_DIFFRACTION_PATTERN
        ///////////////////////////////////
        float cov_model_factor   = cov_model(NdotH, _DiffractionWidth);
        // Cholesky-decomposed cov_init matrix
        float3x3 cov_init_decomp = {_DiffractionCovInit_Row_1.x, _DiffractionCovInit_Row_1.y, _DiffractionCovInit_Row_1.z,
                                    _DiffractionCovInit_Row_2.x, _DiffractionCovInit_Row_2.y, _DiffractionCovInit_Row_2.z,
                                    _DiffractionCovInit_Row_3.x, _DiffractionCovInit_Row_3.y, _DiffractionCovInit_Row_3.z};

        cov_init_decomp *= sqrt(cov_model_factor);
        
        float2 diffraction_uvs   = brdfData.uvs * float2(_DiffractionUV_ScaleX, _DiffractionUV_ScaleY);
#if defined(USING_STEREO_MATRICES)
            float3 cyclopeanEye           = (unity_StereoWorldSpaceCameraPos[0] + unity_StereoWorldSpaceCameraPos[1]) * 0.5;
            float3 cyclopeanViewDirection = SafeNormalize(cyclopeanEye - brdfData.wsPosition);
            float3 cyclopeanHalfDir       = SafeNormalize(lightDirectionWS + cyclopeanViewDirection);
            rnd_number                    = sampleDiffractionPattern(diffraction_uvs, cyclopeanHalfDir, normalWS);
        #else        
            rnd_number                    = sampleDiffractionPattern(diffraction_uvs, halfDir, normalWS);
        #endif

        rnd_number /= 0.15; // In case of OpenSimplex(version 1) implementation, this is the normalization factor. TODO: Is it also the normalization for OpenSimplex2 or Gustavsons simplex ?
        rnd_number  = mul(cov_init_decomp, rnd_number);
        

        // EXPERIMENTAL: Heavy, physically-not-based speckle-approximation
        //const float speckle_cutoff = 0.32;
        //const float speckle_mult   = 10;
        //if (rnd_number.x >= speckle_cutoff)
        //    rnd_number += speckle_mult * (rnd_number.x - speckle_cutoff);
        //////////////////////////////////////////////////////////////////////
#endif

    //float3 F = fresnelSchlickLazanyi(F0, 0, VdotH); // Why should we use Lazanyis improvement to metals ?
    float3 F = fresnelSchlick(F0, 0, VdotH);
    float Ds = distributionGTR2(NdotH, roughness);
    float Gs = geometricSmithGGX(NdotL, roughness) * geometricSmithGGX(NdotV, roughness);
    float3 BRDFspecular = (shift + rnd_number) * (F * Ds * Gs);
    
    return BRDFspecular;
}
#endif




#endif
