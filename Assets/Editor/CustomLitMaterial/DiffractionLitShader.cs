using System;
using System.Collections.Generic;
using UnityEngine;

namespace UnityEditor.Rendering.Universal.ShaderGUI
{
    internal class DiffractionLitShader : BaseShaderGUI
    {
        static readonly string[] workflowModeNames              = Enum.GetNames(typeof(DiffractionLitGUI.WorkflowMode));
        static readonly string[] analyticSpecularBRDFModeNames  = Enum.GetNames(typeof(DiffractionLitGUI.AnalyticSpecularBRDFMode));
        static readonly string[] noiseModeNames                 = Enum.GetNames(typeof(DiffractionLitGUI.NoiseImplementation));

        private DiffractionLitGUI.LitProperties litProperties;
        private DiffractionLitDetailGUI.LitProperties litDetailProperties;

        public override void FillAdditionalFoldouts(MaterialHeaderScopeList materialScopesList)
        {
            materialScopesList.RegisterHeaderScope(DiffractionLitDetailGUI.Styles.detailInputs, Expandable.Details, _ => DiffractionLitDetailGUI.DoDetailArea(litDetailProperties, materialEditor));
        }

        // collect properties from the material properties
        public override void FindProperties(MaterialProperty[] properties)
        {
            base.FindProperties(properties);
            litProperties = new DiffractionLitGUI.LitProperties(properties);
            litDetailProperties = new DiffractionLitDetailGUI.LitProperties(properties);
        }

        // material changed check
        public override void ValidateMaterial(Material material)
        {
            SetMaterialKeywords(material, DiffractionLitGUI.SetMaterialKeywords, DiffractionLitDetailGUI.SetMaterialKeywords);

            if(litProperties.diffractionPattern?.floatValue > 0)
                updateCholeskyDecomposition();
        }

        private void updateCholeskyDecomposition()
        {
            Vector4 column1 = new Vector4(litProperties.diffractionVarR.floatValue,
                                          litProperties.diffractionCovarRG.floatValue,
                                          litProperties.diffractionCovarRB.floatValue,
                                          0);

            Vector4 column2 = new Vector4(litProperties.diffractionCovarRG.floatValue,
                                          litProperties.diffractionVarG.floatValue,
                                          litProperties.diffractionCovarGB.floatValue,
                                          0);

            Vector4 column3 = new Vector4(litProperties.diffractionCovarRB.floatValue,
                                          litProperties.diffractionCovarGB.floatValue,
                                          litProperties.diffractionVarB.floatValue,
                                          0);

            Matrix4x4 cov_init = new Matrix4x4(column1, column2, column3, new Vector4(0,0,0,0));
            Matrix4x4 L        = Matrix4x4.zero;
            

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < (i + 1); j++)
                {
                    float sum = 0;
                    for (int k = 0; k < j; k++)
                        sum += L[i,k] * L[j,k];

                    if (i == j)
                        L[i,j] = Mathf.Sqrt(cov_init[i,i] - sum);
                    else
                        L[i,j] = ((1.0f / L[j,j]) * (cov_init[i,j] - sum));
                }
            }
            /////////////////////
            litProperties.diffractionCovInitMatRow1.vectorValue = new Vector4(L[0, 0], L[0, 1], L[0, 2], 0);
            litProperties.diffractionCovInitMatRow2.vectorValue = new Vector4(L[1, 0], L[1, 1], L[1, 2], 0);
            litProperties.diffractionCovInitMatRow3.vectorValue = new Vector4(L[2, 0], L[2, 1], L[2, 2], 0);
        }

        // material main surface options
        public override void DrawSurfaceOptions(Material material)
        {
            // Use default labelWidth
            EditorGUIUtility.labelWidth = 0f;

            if (litProperties.workflowMode != null)
                DoPopup(DiffractionLitGUI.Styles.workflowModeText, litProperties.workflowMode, workflowModeNames);

            base.DrawSurfaceOptions(material);
        }

        // material main surface inputs
        public override void DrawSurfaceInputs(Material material)
        {
            base.DrawSurfaceInputs(material);
            DiffractionLitGUI.Inputs(litProperties, materialEditor, material);
            DrawEmissionProperties(material, true);
            DrawTileOffset(materialEditor, baseMapProp);
        }

        // material main advanced options
        public override void DrawAdvancedOptions(Material material)
        {
            // Use default labelWidth
            EditorGUIUtility.labelWidth = 0f;

            if (litProperties.analyticSpecularBRDFMode != null)
            {
                DoPopup(DiffractionLitGUI.Styles.analyticSpecularBRDFModeText, litProperties.analyticSpecularBRDFMode, analyticSpecularBRDFModeNames);
                DiffractionLitGUI.AnalyticSpecularBRDFMode specBrdfMode = (DiffractionLitGUI.AnalyticSpecularBRDFMode)litProperties.analyticSpecularBRDFMode.floatValue;
                if (specBrdfMode == DiffractionLitGUI.AnalyticSpecularBRDFMode.CT_DIFF)
                {
                    materialEditor.ShaderProperty(litProperties.diffractionPattern, "Diffraction Pattern");
                    materialEditor.ShaderProperty(litProperties.diffractionWidth,   "Diffraction Width");
                    materialEditor.ShaderProperty(litProperties.diffractionHeight,  "Diffraction Height");

                    materialEditor.ShaderProperty(litProperties.diffractionVarR, "Var_R");
                    materialEditor.ShaderProperty(litProperties.diffractionVarG, "Var_G");
                    materialEditor.ShaderProperty(litProperties.diffractionVarB, "Var_B");
                    materialEditor.ShaderProperty(litProperties.diffractionCovarRG, "Covar_RG");
                    materialEditor.ShaderProperty(litProperties.diffractionCovarRB, "Covar_RB");
                    materialEditor.ShaderProperty(litProperties.diffractionCovarGB, "Covar_GB");

                    materialEditor.ShaderProperty(litProperties.diffractionUVscalingX, "UV Correction X ");
                    materialEditor.ShaderProperty(litProperties.diffractionUVscalingY, "UV Correction Y");
                    materialEditor.ShaderProperty(litProperties.diffractionZWScale, "Spatio-Temporal Pattern-Shift factor");
                    materialEditor.ShaderProperty(litProperties.diffractionStereoSpecularity, "Diffraction Stereoscopic Specularity");
                    DoPopup(DiffractionLitGUI.Styles.noiseModeText, litProperties.noiseImplementation, noiseModeNames);
                }
                    
            }

            if (litProperties.reflections != null && litProperties.highlights != null)
            {
                materialEditor.ShaderProperty(litProperties.highlights, LitGUI.Styles.highlightsText);
                materialEditor.ShaderProperty(litProperties.reflections, LitGUI.Styles.reflectionsText);
            }

            base.DrawAdvancedOptions(material);
        }

        public override void AssignNewShaderToMaterial(Material material, Shader oldShader, Shader newShader)
        {
            if (material == null)
                throw new ArgumentNullException("material");

            // _Emission property is lost after assigning Standard shader to the material
            // thus transfer it before assigning the new shader
            if (material.HasProperty("_Emission"))
            {
                material.SetColor("_EmissionColor", material.GetColor("_Emission"));
            }

            base.AssignNewShaderToMaterial(material, oldShader, newShader);

            if (oldShader == null || !oldShader.name.Contains("Legacy Shaders/"))
            {
                SetupMaterialBlendMode(material);
                return;
            }

            SurfaceType surfaceType = SurfaceType.Opaque;
            BlendMode blendMode = BlendMode.Alpha;
            if (oldShader.name.Contains("/Transparent/Cutout/"))
            {
                surfaceType = SurfaceType.Opaque;
                material.SetFloat("_AlphaClip", 1);
            }
            else if (oldShader.name.Contains("/Transparent/"))
            {
                // NOTE: legacy shaders did not provide physically based transparency
                // therefore Fade mode
                surfaceType = SurfaceType.Transparent;
                blendMode = BlendMode.Alpha;
            }
            material.SetFloat("_Blend", (float)blendMode);

            material.SetFloat("_Surface", (float)surfaceType);
            if (surfaceType == SurfaceType.Opaque)
            {
                material.DisableKeyword("_SURFACE_TYPE_TRANSPARENT");
            }
            else
            {
                material.EnableKeyword("_SURFACE_TYPE_TRANSPARENT");
            }

            if (oldShader.name.Equals("Standard (Specular setup)"))
            {
                material.SetFloat("_WorkflowMode", (float)DiffractionLitGUI.WorkflowMode.Specular);
                Texture texture = material.GetTexture("_SpecGlossMap");
                if (texture != null)
                    material.SetTexture("_MetallicSpecGlossMap", texture);
            }
            else
            {
                material.SetFloat("_WorkflowMode", (float)DiffractionLitGUI.WorkflowMode.Metallic);
                Texture texture = material.GetTexture("_MetallicGlossMap");
                if (texture != null)
                    material.SetTexture("_MetallicSpecGlossMap", texture);
            }
        }
    }
}
