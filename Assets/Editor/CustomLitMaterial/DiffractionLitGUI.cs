using UnityEngine;
using UnityEngine.Rendering;

namespace UnityEditor.Rendering.Universal.ShaderGUI
{

    internal static class Property
    {
        public static readonly string SpecularWorkflowMode = "_WorkflowMode";
        public static readonly string SurfaceType = "_Surface";
        public static readonly string BlendMode = "_Blend";
        public static readonly string AlphaClip = "_AlphaClip";
        public static readonly string SrcBlend = "_SrcBlend";
        public static readonly string DstBlend = "_DstBlend";
        public static readonly string ZWrite = "_ZWrite";
        public static readonly string CullMode = "_Cull";
        public static readonly string CastShadows = "_CastShadows";
        public static readonly string ReceiveShadows = "_ReceiveShadows";
        public static readonly string QueueOffset = "_QueueOffset";

        // for ShaderGraph shaders only
        public static readonly string ZTest = "_ZTest";
        public static readonly string ZWriteControl = "_ZWriteControl";
        public static readonly string QueueControl = "_QueueControl";

        // Global Illumination requires some properties to be named specifically:
        public static readonly string EmissionMap = "_EmissionMap";
        public static readonly string EmissionColor = "_EmissionColor";
    }

    public static class DiffractionLitGUI
    {
        public enum WorkflowMode
        {
            Specular = 0,
            Metallic
        }

        public enum AnalyticSpecularBRDFMode
        {
            Unity_Default = 0, 
            CT = 1, 
            CT_DIFF = 2
        }

        public enum NoiseImplementation
        {
            GUSTAVSON_SIMPLEX = 0,
            OPEN_SIMPLEX_2 = 1
        }

        public enum SmoothnessMapChannel
        {
            SpecularMetallicAlpha,
            AlbedoAlpha,
        }

        public static class Styles
        {

	    public static GUIContent analyticSpecularBRDFModeText = EditorGUIUtility.TrTextContent("Analytic Specular BRDF",
		"Select the which BRDF model is used for the rendering of specular reflections from analytic(direct) light sources.");

            public static GUIContent noiseTextureDiffractionText = EditorGUIUtility.TrTextContent("Diffraction Distribution",
                "Noise distribution used for the diffraction pattern.");

            public static GUIContent noiseModeText = EditorGUIUtility.TrTextContent("Diffraction Noise",
                "Select between the simplex noise implementation of Steven Gustavson or Open Simplex 2");

            public static GUIContent workflowModeText = EditorGUIUtility.TrTextContent("Workflow Mode",
                "Select a workflow that fits your textures. Choose between Metallic or Specular.");

            public static GUIContent specularMapText =
                EditorGUIUtility.TrTextContent("Specular Map", "Designates a Specular Map and specular color determining the apperance of reflections on this Material's surface.");

            public static GUIContent metallicMapText =
                EditorGUIUtility.TrTextContent("Metallic Map", "Sets and configures the map for the Metallic workflow.");

            public static GUIContent smoothnessText = EditorGUIUtility.TrTextContent("Smoothness",
                "Controls the spread of highlights and reflections on the surface.");

            public static GUIContent smoothnessMapChannelText =
                EditorGUIUtility.TrTextContent("Source",
                    "Specifies where to sample a smoothness map from. By default, uses the alpha channel for your map.");

            public static GUIContent highlightsText = EditorGUIUtility.TrTextContent("Specular Highlights",
                "When enabled, the Material reflects the shine from direct lighting.");

            public static GUIContent reflectionsText =
                EditorGUIUtility.TrTextContent("Environment Reflections",
                    "When enabled, the Material samples reflections from the nearest Reflection Probes or Lighting Probe.");
  
            public static GUIContent heightMapText = EditorGUIUtility.TrTextContent("Height Map",
                "Defines a Height Map that will drive a parallax effect in the shader making the surface seem displaced.");

            public static GUIContent occlusionText = EditorGUIUtility.TrTextContent("Occlusion Map",
                "Sets an occlusion map to simulate shadowing from ambient lighting.");

            public static readonly string[] metallicSmoothnessChannelNames = { "Metallic Alpha", "Albedo Alpha" };
            public static readonly string[] specularSmoothnessChannelNames = { "Specular Alpha", "Albedo Alpha" };

            public static GUIContent clearCoatText = EditorGUIUtility.TrTextContent("Clear Coat",
                "A multi-layer material feature which simulates a thin layer of coating on top of the surface material." +
                "\nPerformance cost is considerable as the specular component is evaluated twice, once per layer.");

            public static GUIContent clearCoatMaskText = EditorGUIUtility.TrTextContent("Mask",
                "Specifies the amount of the coat blending." +
                "\nActs as a multiplier of the clear coat map mask value or as a direct mask value if no map is specified." +
                "\nThe map specifies clear coat mask in the red channel and clear coat smoothness in the green channel.");

            public static GUIContent clearCoatSmoothnessText = EditorGUIUtility.TrTextContent("Smoothness",
                "Specifies the smoothness of the coating." +
                "\nActs as a multiplier of the clear coat map smoothness value or as a direct smoothness value if no map is specified.");
        }

        public struct LitProperties
        {
            // Surface Option Props
            public MaterialProperty workflowMode;

            // Surface Input Props
            public MaterialProperty metallic;
            public MaterialProperty specColor;
            public MaterialProperty metallicGlossMap;
            public MaterialProperty specGlossMap;
            public MaterialProperty smoothness;
            public MaterialProperty smoothnessMapChannel;
            public MaterialProperty bumpMapProp;
            public MaterialProperty bumpScaleProp;
            public MaterialProperty parallaxMapProp;
            public MaterialProperty parallaxScaleProp;
            public MaterialProperty occlusionStrength;
            public MaterialProperty occlusionMap;

            // Advanced Props
            public MaterialProperty highlights;
            public MaterialProperty reflections;
            public MaterialProperty analyticSpecularBRDFMode;
            public MaterialProperty noiseImplementation;
            public MaterialProperty diffractionPattern;
            public MaterialProperty diffractionWidth;
            public MaterialProperty diffractionHeight;

            public MaterialProperty diffractionVarR;
            public MaterialProperty diffractionVarG;
            public MaterialProperty diffractionVarB;
            public MaterialProperty diffractionCovarRG;
            public MaterialProperty diffractionCovarRB;
            public MaterialProperty diffractionCovarGB;

            public MaterialProperty diffractionCovInitMatRow1;
            public MaterialProperty diffractionCovInitMatRow2;
            public MaterialProperty diffractionCovInitMatRow3;

            public MaterialProperty diffractionStereoSpecularity;
            public MaterialProperty diffractionZWScale;
            public MaterialProperty diffractionUVscalingX;
            public MaterialProperty diffractionUVscalingY;
            

            public MaterialProperty clearCoat;  // Enable/Disable dummy property
            public MaterialProperty clearCoatMap;
            public MaterialProperty clearCoatMask;
            public MaterialProperty clearCoatSmoothness;

            public LitProperties(MaterialProperty[] properties)
            {
                // Surface Option Props
                workflowMode = BaseShaderGUI.FindProperty("_WorkflowMode", properties, false);
                // Surface Input Props
                metallic = BaseShaderGUI.FindProperty("_Metallic", properties);
                specColor = BaseShaderGUI.FindProperty("_SpecColor", properties, false);
                metallicGlossMap = BaseShaderGUI.FindProperty("_MetallicGlossMap", properties);
                specGlossMap = BaseShaderGUI.FindProperty("_SpecGlossMap", properties, false);
                smoothness = BaseShaderGUI.FindProperty("_Smoothness", properties, false);
                smoothnessMapChannel = BaseShaderGUI.FindProperty("_SmoothnessTextureChannel", properties, false);
                bumpMapProp = BaseShaderGUI.FindProperty("_BumpMap", properties, false);
                bumpScaleProp = BaseShaderGUI.FindProperty("_BumpScale", properties, false);
                parallaxMapProp = BaseShaderGUI.FindProperty("_ParallaxMap", properties, false);
                parallaxScaleProp = BaseShaderGUI.FindProperty("_Parallax", properties, false);
                occlusionStrength = BaseShaderGUI.FindProperty("_OcclusionStrength", properties, false);
                occlusionMap = BaseShaderGUI.FindProperty("_OcclusionMap", properties, false);
                // Advanced Props
                highlights = BaseShaderGUI.FindProperty("_SpecularHighlights", properties, false);
                reflections = BaseShaderGUI.FindProperty("_EnvironmentReflections", properties, false);
		        analyticSpecularBRDFMode     = BaseShaderGUI.FindProperty("_AnalyticSpecularBRDFMode", properties, false);
                noiseImplementation          = BaseShaderGUI.FindProperty("_NoiseImplementationDiffraction", properties, false);
                diffractionWidth             = BaseShaderGUI.FindProperty("_DiffractionWidth", properties, false);
                diffractionHeight            = BaseShaderGUI.FindProperty("_DiffractionHeight", properties, false);

                diffractionVarR              = BaseShaderGUI.FindProperty("_DiffractionVar_R", properties, false);
                diffractionVarG              = BaseShaderGUI.FindProperty("_DiffractionVar_G", properties, false);
                diffractionVarB              = BaseShaderGUI.FindProperty("_DiffractionVar_B", properties, false);
                diffractionCovarRG           = BaseShaderGUI.FindProperty("_DiffractionCovar_RG", properties, false);
                diffractionCovarRB           = BaseShaderGUI.FindProperty("_DiffractionCovar_RB", properties, false);
                diffractionCovarGB           = BaseShaderGUI.FindProperty("_DiffractionCovar_GB", properties, false);
                diffractionCovInitMatRow1    = BaseShaderGUI.FindProperty("_DiffractionCovInit_Row_1", properties, false);
                diffractionCovInitMatRow2    = BaseShaderGUI.FindProperty("_DiffractionCovInit_Row_2", properties, false);
                diffractionCovInitMatRow3    = BaseShaderGUI.FindProperty("_DiffractionCovInit_Row_3", properties, false);

                diffractionPattern           = BaseShaderGUI.FindProperty("_DiffractionPatternToggle", properties, false);
                diffractionStereoSpecularity = BaseShaderGUI.FindProperty("_DiffractionStereoSpecularity", properties, false);
                diffractionZWScale           = BaseShaderGUI.FindProperty("_DiffractionZW_Scale", properties, false);
                diffractionUVscalingX        = BaseShaderGUI.FindProperty("_DiffractionUV_ScaleX", properties, false);
                diffractionUVscalingY        = BaseShaderGUI.FindProperty("_DiffractionUV_ScaleY", properties, false);

                clearCoat = BaseShaderGUI.FindProperty("_ClearCoat", properties, false);
                clearCoatMap = BaseShaderGUI.FindProperty("_ClearCoatMap", properties, false);
                clearCoatMask = BaseShaderGUI.FindProperty("_ClearCoatMask", properties, false);
                clearCoatSmoothness = BaseShaderGUI.FindProperty("_ClearCoatSmoothness", properties, false);
            }
        }

        public static void Inputs(LitProperties properties, MaterialEditor materialEditor, Material material)
        {
            DoMetallicSpecularArea(properties, materialEditor, material);
            BaseShaderGUI.DrawNormalArea(materialEditor, properties.bumpMapProp, properties.bumpScaleProp);

            if (HeightmapAvailable(material))
                DoHeightmapArea(properties, materialEditor);

            if (properties.occlusionMap != null)
            {
                materialEditor.TexturePropertySingleLine(Styles.occlusionText, properties.occlusionMap,
                    properties.occlusionMap.textureValue != null ? properties.occlusionStrength : null);
            }

            // Check that we have all the required properties for clear coat,
            // otherwise we will get null ref exception from MaterialEditor GUI helpers.
            if (ClearCoatAvailable(material))
                DoClearCoat(properties, materialEditor, material);
        }

        private static bool ClearCoatAvailable(Material material)
        {
            return material.HasProperty("_ClearCoat")
                && material.HasProperty("_ClearCoatMap")
                && material.HasProperty("_ClearCoatMask")
                && material.HasProperty("_ClearCoatSmoothness");
        }

        private static bool HeightmapAvailable(Material material)
        {
            return material.HasProperty("_Parallax")
                && material.HasProperty("_ParallaxMap");
        }

        private static void DoHeightmapArea(LitProperties properties, MaterialEditor materialEditor)
        {
            materialEditor.TexturePropertySingleLine(Styles.heightMapText, properties.parallaxMapProp,
                properties.parallaxMapProp.textureValue != null ? properties.parallaxScaleProp : null);
        }

        private static bool ClearCoatEnabled(Material material)
        {
            return material.HasProperty("_ClearCoat") && material.GetFloat("_ClearCoat") > 0.0;
        }

        public static void DoClearCoat(LitProperties properties, MaterialEditor materialEditor, Material material)
        {
            materialEditor.ShaderProperty(properties.clearCoat, Styles.clearCoatText);
            var coatEnabled = material.GetFloat("_ClearCoat") > 0.0;

            EditorGUI.BeginDisabledGroup(!coatEnabled);
            {
                materialEditor.TexturePropertySingleLine(Styles.clearCoatMaskText, properties.clearCoatMap, properties.clearCoatMask);

                EditorGUI.indentLevel += 2;

                // Texture and HDR color controls
                materialEditor.ShaderProperty(properties.clearCoatSmoothness, Styles.clearCoatSmoothnessText);

                EditorGUI.indentLevel -= 2;
            }
            EditorGUI.EndDisabledGroup();
        }

        public static void DoMetallicSpecularArea(LitProperties properties, MaterialEditor materialEditor, Material material)
        {
            string[] smoothnessChannelNames;
            bool hasGlossMap = false;
            if (properties.workflowMode == null ||
                (WorkflowMode)properties.workflowMode.floatValue == WorkflowMode.Metallic)
            {
                hasGlossMap = properties.metallicGlossMap.textureValue != null;
                smoothnessChannelNames = Styles.metallicSmoothnessChannelNames;
                materialEditor.TexturePropertySingleLine(Styles.metallicMapText, properties.metallicGlossMap,
                    hasGlossMap ? null : properties.metallic);
            }
            else
            {
                hasGlossMap = properties.specGlossMap.textureValue != null;
                smoothnessChannelNames = Styles.specularSmoothnessChannelNames;
                BaseShaderGUI.TextureColorProps(materialEditor, Styles.specularMapText, properties.specGlossMap,
                    hasGlossMap ? null : properties.specColor);
            }
            DoSmoothness(materialEditor, material, properties.smoothness, properties.smoothnessMapChannel, smoothnessChannelNames);
        }

        internal static bool IsOpaque(Material material)
        {
            bool opaque = true;
            if (material.HasProperty(Property.SurfaceType))
                opaque = ((BaseShaderGUI.SurfaceType)material.GetFloat(Property.SurfaceType) == BaseShaderGUI.SurfaceType.Opaque);
            return opaque;
        }

        public static void DoSmoothness(MaterialEditor materialEditor, Material material, MaterialProperty smoothness, MaterialProperty smoothnessMapChannel, string[] smoothnessChannelNames)
        {
            EditorGUI.indentLevel += 2;

            materialEditor.ShaderProperty(smoothness, Styles.smoothnessText);

            if (smoothnessMapChannel != null) // smoothness channel
            {
                var opaque = IsOpaque(material);
                EditorGUI.indentLevel++;
                EditorGUI.showMixedValue = smoothnessMapChannel.hasMixedValue;
                if (opaque)
                {
                    EditorGUI.BeginChangeCheck();
                    var smoothnessSource = (int)smoothnessMapChannel.floatValue;
                    smoothnessSource = EditorGUILayout.Popup(Styles.smoothnessMapChannelText, smoothnessSource, smoothnessChannelNames);
                    if (EditorGUI.EndChangeCheck())
                        smoothnessMapChannel.floatValue = smoothnessSource;
                }
                else
                {
                    EditorGUI.BeginDisabledGroup(true);
                    EditorGUILayout.Popup(Styles.smoothnessMapChannelText, 0, smoothnessChannelNames);
                    EditorGUI.EndDisabledGroup();
                }
                EditorGUI.showMixedValue = false;
                EditorGUI.indentLevel--;
            }
            EditorGUI.indentLevel -= 2;
        }

        public static SmoothnessMapChannel GetSmoothnessMapChannel(Material material)
        {
            int ch = (int)material.GetFloat("_SmoothnessTextureChannel");
            if (ch == (int)SmoothnessMapChannel.AlbedoAlpha)
                return SmoothnessMapChannel.AlbedoAlpha;

            return SmoothnessMapChannel.SpecularMetallicAlpha;
        }

        // (shared by all lit shaders, including shadergraph Lit Target and Lit.shader)
        internal static void SetupSpecularWorkflowKeyword(Material material, out bool isSpecularWorkflow)
        {
            isSpecularWorkflow = false;     // default is metallic workflow
            if (material.HasProperty(Property.SpecularWorkflowMode))
                isSpecularWorkflow = ((WorkflowMode)material.GetFloat(Property.SpecularWorkflowMode)) == WorkflowMode.Specular;
            CoreUtils.SetKeyword(material, "_SPECULAR_SETUP", isSpecularWorkflow);
        }

        // setup keywords for Lit.shader
        public static void SetMaterialKeywords(Material material)
        {
            SetupSpecularWorkflowKeyword(material, out bool isSpecularWorkFlow);

            // Note: keywords must be based on Material value not on MaterialProperty due to multi-edit & material animation
            // (MaterialProperty value might come from renderer material property block)
            var specularGlossMap = isSpecularWorkFlow ? "_SpecGlossMap" : "_MetallicGlossMap";
            var hasGlossMap = material.GetTexture(specularGlossMap) != null;

            CoreUtils.SetKeyword(material, "_METALLICSPECGLOSSMAP", hasGlossMap);

            if (material.HasProperty("_SpecularHighlights"))
                CoreUtils.SetKeyword(material, "_SPECULARHIGHLIGHTS_OFF",
                    material.GetFloat("_SpecularHighlights") == 0.0f);
            if (material.HasProperty("_EnvironmentReflections"))
                CoreUtils.SetKeyword(material, "_ENVIRONMENTREFLECTIONS_OFF",
                    material.GetFloat("_EnvironmentReflections") == 0.0f);
            if (material.HasProperty("_OcclusionMap"))
                CoreUtils.SetKeyword(material, "_OCCLUSIONMAP", material.GetTexture("_OcclusionMap"));

            if (material.HasProperty("_ParallaxMap"))
                CoreUtils.SetKeyword(material, "_PARALLAXMAP", material.GetTexture("_ParallaxMap"));

            if (material.HasProperty("_SmoothnessTextureChannel"))
            {
                var opaque = IsOpaque(material);
                CoreUtils.SetKeyword(material, "_SMOOTHNESS_TEXTURE_ALBEDO_CHANNEL_A",
                    GetSmoothnessMapChannel(material) == SmoothnessMapChannel.AlbedoAlpha && opaque);
            }

            if (material.HasProperty("_AnalyticSpecularBRDFMode"))
            {
                AnalyticSpecularBRDFMode specBrdfMode = (AnalyticSpecularBRDFMode)material.GetFloat("_AnalyticSpecularBRDFMode");
                switch (specBrdfMode)
                {
                    case AnalyticSpecularBRDFMode.Unity_Default:
                        CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_COOK_TORRANCE", false);
                        CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_DIFFRACTION", false);
                        break;
                    case AnalyticSpecularBRDFMode.CT:
                        CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_COOK_TORRANCE", true);
                        CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_DIFFRACTION", false);
                        break;
                    case AnalyticSpecularBRDFMode.CT_DIFF:
                        CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_COOK_TORRANCE", false);
                        CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_DIFFRACTION", true);
                        break;
                }

                if (material.HasProperty("_DiffractionPatternToggle"))
                {
                    CoreUtils.SetKeyword(material, "_DIRECT_LIGHT_BRDF_DIFFRACTION_PATTERN",
                    material.GetFloat("_DiffractionPatternToggle") == 1.0f);
                }

                if (material.HasProperty("_NoiseImplementationDiffraction"))
                    CoreUtils.SetKeyword(material, "_DIFFRACTION_PATTERN_OPEN_SIMPLEX_2",
                    material.GetFloat("_NoiseImplementationDiffraction") == 1.0f);
            }

            // Clear coat keywords are independent to remove possiblity of invalid combinations.
            if (ClearCoatEnabled(material))
            {
                var hasMap = material.HasProperty("_ClearCoatMap") && material.GetTexture("_ClearCoatMap") != null;
                if (hasMap)
                {
                    CoreUtils.SetKeyword(material, "_CLEARCOAT", false);
                    CoreUtils.SetKeyword(material, "_CLEARCOATMAP", true);
                }
                else
                {
                    CoreUtils.SetKeyword(material, "_CLEARCOAT", true);
                    CoreUtils.SetKeyword(material, "_CLEARCOATMAP", false);
                }
            }
            else
            {
                CoreUtils.SetKeyword(material, "_CLEARCOAT", false);
                CoreUtils.SetKeyword(material, "_CLEARCOATMAP", false);
            }
        }
    }
}
