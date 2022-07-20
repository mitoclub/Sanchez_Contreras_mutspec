import sys

rule all:
    input:
        "figures/Figure_1.png",
        "figures/Figure_2.png",
        "figures/Figure_3.png",
        "figures/Figure_4.png",
        "figures/Figure_5.png",
        "figures/Figure_6.png",
        "figures/Supp_Figure_1.png",
        "figures/Supp_Figure_2.png",
        "figures/Supp_Figure_3.png",
        "figures/Supp_Figure_4.png",
        "figures/Supp_Figure_8.png",
        "figures/Supp_Figure_9.png",
        "figures/Supp_Figure_10.png",
        "figures/Supp_Figure_11.png",
    
        "data/stats/Figure_1A_statistics.csv",
        "data/stats/Figure_1D_statistics.csv",
        "data/stats/Figure_1A_Old_heatmap_stats.csv",
        "data/stats/Figure_1A_Young_heatmap_stats.csv",
        "data/stats/Figure_1D_Old_heatmap_stats.csv",
        "data/stats/Figure_1D_Young_heatmap_stats.csv",
        "data/stats/Figure_2_Old_AC_TG_heatmap_stats.csv",
        "data/stats/Figure_2_Old_AG_TC_heatmap_stats.csv",
        "data/stats/Figure_2_Old_AT_TA_heatmap_stats.csv",
        "data/stats/Figure_2_Old_CA_GT_heatmap_stats.csv",
        "data/stats/Figure_2_Old_CG_GC_heatmap_stats.csv",
        "data/stats/Figure_2_Old_CT_GA_heatmap_stats.csv",
        "data/stats/Figure_2_Young_AC_TG_heatmap_stats.csv",
        "data/stats/Figure_2_Young_AG_TC_heatmap_stats.csv",
        "data/stats/Figure_2_Young_AT_TA_heatmap_stats.csv",
        "data/stats/Figure_2_Young_CA_GT_heatmap_stats.csv",
        "data/stats/Figure_2_Young_CG_GC_heatmap_stats.csv",
        "data/stats/Figure_2_Young_CT_GA_heatmap_stats.csv",
        "data/stats/Figure_3_ratio_statistics.csv",
        "data/stats/Figure_4A_statistics.csv",
        "data/stats/Figure_4B_statistics.csv",
        "data/stats/Figure_6_Dunnett_statistics.csv",
        "data/stats/Supplemental_Figure_8_statistics.csv",
        
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv",
        "data/imported_data/mut_file_data.csv",
        "data/imported_data/summary_clone_data.csv",
    output:
        temp(touch(".ruleAllFinished"))

rule initializeEnvs:
    input:
        "full.initialized"
    params:
        basePath = sys.path[0]
    output:
        temp(".env_initialized")
    shell:
        """
        touch "{output}"
        """

rule initializeFullEnv:
    output:
        temp("full.initialized")
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        echo "Initializing envs"
        touch full.initialized
        """
rule makeDirs:
    output:
        temp(touch(".dirsMade")) 
    shell:
        """
        mkdir -p data/stats \
        data/imported_data \
        figures
        """
rule compile_data:
    output:
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv",
        "data/imported_data/mut_file_data.csv",
        "data/imported_data/summary_clone_data.csv"
    priority:
        1
    shell:
        """
        python compile_data.py
        """

rule compute_stats:
    input:
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv",
        "data/imported_data/mut_file_data.csv",
        "data/imported_data/summary_clone_data.csv",
    output:
        "data/stats/Figure_1A_statistics.csv",
        "data/stats/Figure_1D_statistics.csv",
        "data/stats/Figure_1A_Old_heatmap_stats.csv",
        "data/stats/Figure_1A_Young_heatmap_stats.csv",
        "data/stats/Figure_1D_Old_heatmap_stats.csv",
        "data/stats/Figure_1D_Young_heatmap_stats.csv",
        "data/stats/Figure_2_Old_AC_TG_heatmap_stats.csv",
        "data/stats/Figure_2_Old_AG_TC_heatmap_stats.csv",
        "data/stats/Figure_2_Old_AT_TA_heatmap_stats.csv",
        "data/stats/Figure_2_Old_CA_GT_heatmap_stats.csv",
        "data/stats/Figure_2_Old_CG_GC_heatmap_stats.csv",
        "data/stats/Figure_2_Old_CT_GA_heatmap_stats.csv",
        "data/stats/Figure_2_Young_AC_TG_heatmap_stats.csv",
        "data/stats/Figure_2_Young_AG_TC_heatmap_stats.csv",
        "data/stats/Figure_2_Young_AT_TA_heatmap_stats.csv",
        "data/stats/Figure_2_Young_CA_GT_heatmap_stats.csv",
        "data/stats/Figure_2_Young_CG_GC_heatmap_stats.csv",
        "data/stats/Figure_2_Young_CT_GA_heatmap_stats.csv",
        "data/stats/Figure_3_ratio_statistics.csv",
        "data/stats/Figure_4A_statistics.csv",
        "data/stats/Figure_4B_statistics.csv",
        "data/stats/Figure_6_Dunnett_statistics.csv",
        "data/stats/Supplemental_Figure_8_statistics.csv"
    priority:
        1
    shell:
        """
        python compute_stats.py
        """
        
rule Figure_1:
    input:
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv"
    output: 
        "figures/Figure_1.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Figure_1.py
        """
rule Figure_2:
    input:
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv"
    output:
        "figures/Figure_2.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Figure_2.py
        """

rule Figure_3:
    input:
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv"
    output:
        "figures/Figure_3.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Figure_3.py
        """

rule Figure_4:
    input:
        "data/imported_data/mut_file_data.csv",
        "data/imported_data/summary_clone_data.csv"
    output:
        "figures/Figure_4.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Figure_4.py
        """

rule Figure_5:
    input:
        "data/imported_data/mut_file_data.csv",
        "data/imported_data/summary_clone_data.csv"
    output:
        "figures/Figure_5.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Figure_5.py
        """

rule Figure_6:
    input:
        "data/imported_data/summary_data_wide.csv",
        "data/imported_data/summary_data_tidy.csv"
    output:
        "figures/Figure_6.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Figure_6.py
        """

rule Supplemental_Figure_1:
    output:
        "figures/Supp_Figure_1.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_1.py
        """

rule Supplemental_Figure_2:
    output:
        "figures/Supp_Figure_2.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_2.py
        """
        
rule Supplemental_Figure_3:
    output:
        "figures/Supp_Figure_3.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_3.py
        """
        
rule Supplemental_Figure_4:
    output:
        "figures/Supp_Figure_4.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_4.py
        """
        
rule Supplemental_Figure_8:
    output:
        "figures/Supp_Figure_8.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_8.py
        """
        
rule Supplemental_Figure_9:
    output:
        "figures/Supp_Figure_9.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_9.py
        """
        
rule Supplemental_Figure_10:
    output:
        "figures/Supp_Figure_10.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_10.py
        """
        
rule Supplemental_Figure_11:
    output:
        "figures/Supp_Figure_11.png"
    conda:
        "Mouse_mtDNA_analysis_env.yaml"
    shell:
        """
        python3 Supp_Figure_11.py
        """