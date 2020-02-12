#!/usr/bin/env nextflow

import groovy.json.*

params.root = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["b0_thr_extract_b0":"$params.b0_thr_extract_b0",
                "dwi_shell_tolerance":"$params.dwi_shell_tolerance",
                "dilate_b0_mask_prelim_brain_extraction":"$params.dilate_b0_mask_prelim_brain_extraction",
                "bet_prelim_f":"$params.bet_prelim_f",
                "run_dwi_denoising":"$params.run_dwi_denoising",
                "extent":"$params.extent",
                "run_topup":"$params.run_topup",
                "encoding_direction":"$params.encoding_direction",
                "readout":"$params.readout",
                "run_eddy":"$params.run_eddy",
                "eddy_cmd":"$params.eddy_cmd",
                "bet_topup_before_eddy_f":"$params.bet_topup_before_eddy_f",
                "use_slice_drop_correction":"$params.use_slice_drop_correction",
                "bet_dwi_final_f":"$params.bet_dwi_final_f",
                "run_resample_dwi":"$params.run_resample_dwi",
                "dwi_resolution":"$params.dwi_resolution",
                "dwi_interpolation":"$params.dwi_interpolation",
                "cpu_count":"$cpu_count",
                "processes_denoise_dwi":"$params.processes_denoise_dwi",
                "processes_eddy":"$params.processes_eddy"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Preprocessing pipeline"
log.info "======================"
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Denoise DWI]"
log.info "Denoise DWI: $params.run_dwi_denoising"
log.info ""
log.info "[Topup]"
log.info "Run Topup: $params.run_topup"
log.info ""
log.info "[Eddy]"
log.info "Run Eddy: $params.run_eddy"
log.info "Eddy command: $params.eddy_cmd"
log.info ""
log.info "[Resample DWI]"
log.info "Resample DWI: $params.run_resample_dwi"
log.info "Resolution: $params.dwi_resolution"

log.info "Number of processes per tasks"
log.info "============================="
log.info "Denoise DWI: $params.processes_denoise_dwi"
log.info "Eddy: $params.processes_eddy"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    data = Channel
        .fromFilePairs("$root/**/*{bval,bvec,dwi.nii.gz}",
                       size: 3,
                       maxDepth:1,
                       flat: true) {it.parent.name}

    data
        .map{[it, params.readout, params.encoding_direction].flatten()}
        .into{in_data; check_subjects_number}

    Channel
    .fromPath("$root/**/*rev_b0.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{rev_b0; check_rev_b0}
}
else {
    error "Error ~ Please use --root, --bids or --bids_config for the input data."
}

(dwi, gradients, readout_encoding) = in_data
    .map{sid, bvals, bvecs, dwi, readout, encoding -> [tuple(sid, dwi),
                                        tuple(sid, bvals, bvecs),
                                        tuple(sid, readout, encoding)]}
    .separate(3)

check_rev_b0.count().into{ rev_b0_counter; number_rev_b0_for_compare }

check_subjects_number.count().into{ number_subj_for_null_check; number_subj_for_compare }

if (number_subj_for_null_check == 0){
    error "Error ~ No subjects found. Please check the naming convention or your BIDS folder."
}

number_subj_for_compare.count()
    .concat(number_rev_b0_for_compare)
    .toList()
    .subscribe{a, b -> if (a != b && b > 0) 
    error "Error ~ Some subjects have a reversed phase encoded b=0 and others don't.\n" +
          "Please be sure to have the same acquisitions for all subjects."}

dwi.into{dwi_for_prelim_bet; dwi_for_denoise}

gradients
    .into{gradients_for_prelim_bet; gradients_for_eddy; gradients_for_topup;
          gradients_for_eddy_topup}

readout_encoding
    .into{readout_encoding_for_topup; readout_encoding_for_eddy;
          readout_encoding_for_eddy_topup}

dwi_for_prelim_bet
    .join(gradients_for_prelim_bet)
    .set{dwi_gradient_for_prelim_bet}

process README {
    cpus 1
    publishDir = params.Readme_Publish_Dir
    tag = "README"

    output:
    file "readme.txt"

    script:
    String list_options = new String();
    for (String item : params) {
        list_options += item + "\n"
    }
    """
    echo "TractoFlow pipeline\n" >> readme.txt
    echo "Start time: $workflow.start\n" >> readme.txt
    echo "[Command-line]\n$workflow.commandLine\n" >> readme.txt
    echo "[Git Info]\n" >> readme.txt
    echo "$workflow.repository - $workflow.revision [$workflow.commitId]\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "$list_options" >> readme.txt
    """
}

process Bet_Prelim_DWI {
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradient_for_prelim_bet
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__b0_bet_mask_dilated.nii.gz" into\
        b0_mask_for_eddy
    file "${sid}__b0_bet.nii.gz"
    file "${sid}__b0_bet_mask.nii.gz"

    when:
    rev_b0_count == 0 || !params.run_topup || (!params.run_eddy && params.run_topup)

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    bet ${sid}__b0.nii.gz ${sid}__b0_bet.nii.gz -m -R -f $params.bet_prelim_f
    maskfilter ${sid}__b0_bet_mask.nii.gz dilate ${sid}__b0_bet_mask_dilated.nii.gz\
        --npass $params.dilate_b0_mask_prelim_brain_extraction -nthreads 1
    mrcalc ${sid}__b0.nii.gz ${sid}__b0_bet_mask_dilated.nii.gz\
        -mult ${sid}__b0_bet.nii.gz -quiet -force -nthreads 1
    """
}

process Denoise_DWI {
    cpus params.processes_denoise_dwi

    input:
    set sid, file(dwi) from dwi_for_denoise

    output:
    set sid, "${sid}__dwi_denoised.nii.gz" into\
        dwi_for_eddy,
        dwi_for_topup,
        dwi_for_eddy_topup

    script:
    // The denoised DWI is clipped to 0 since negative values
    // could have been introduced.
    if(params.run_dwi_denoising)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        dwidenoise $dwi dwi_denoised.nii.gz -extent $params.extent -nthreads $task.cpus
        fslmaths dwi_denoised.nii.gz -thr 0 ${sid}__dwi_denoised.nii.gz
        """
    else
        """
        mv $dwi ${sid}__dwi_denoised.nii.gz
        """
}

dwi_for_topup
    .join(gradients_for_topup)
    .join(rev_b0)
    .join(readout_encoding_for_topup)
    .set{dwi_gradients_rev_b0_for_topup}

process Topup {
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec), file(rev_b0), readout, encoding\
        from dwi_gradients_rev_b0_for_topup

    output:
    set sid, "${sid}__corrected_b0s.nii.gz", "${params.prefix_topup}_fieldcoef.nii.gz",
    "${params.prefix_topup}_movpar.txt" into topup_files_for_eddy_topup
    file "${sid}__rev_b0_warped.nii.gz"

    when:
    params.run_topup && params.run_eddy

    script:
    """
    export OMP_NUM_THREADS=$task.cpus
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    extract_b0.py $dwi $bval $bvec b0_mean.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    antsRegistrationSyNQuick.sh -d 3 -f b0_mean.nii.gz -m $rev_b0 -o output -t r -e 1
    mv outputWarped.nii.gz ${sid}__rev_b0_warped.nii.gz
    scil_prepare_topup_command.py $dwi $bval $bvec ${sid}__rev_b0_warped.nii.gz\
        --config $params.config_topup --b0_thr $params.b0_thr_extract_b0\
        --encoding_direction $encoding\
        --readout $readout --output_prefix $params.prefix_topup\
        --output_script
    sh topup.sh
    cp corrected_b0s.nii.gz ${sid}__corrected_b0s.nii.gz
    """
}

dwi_for_eddy
    .join(gradients_for_eddy)
    .join(b0_mask_for_eddy)
    .join(readout_encoding_for_eddy)
    .set{dwi_gradients_mask_topup_files_for_eddy}

process Eddy {
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), file(mask), readout, encoding\
        from dwi_gradients_mask_topup_files_for_eddy
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__dwi_corrected.nii.gz", "${sid}__bval_eddy",
        "${sid}__dwi_eddy_corrected.bvec" into\
        dwi_gradients_from_eddy
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_from_eddy
    set sid, "${sid}__bval_eddy", "${sid}__dwi_eddy_corrected.bvec" into\
        gradients_from_eddy

    when:
    rev_b0_count == 0 || !params.run_topup || (!params.run_eddy && params.run_topup)

    // Corrected DWI is clipped to 0 since Eddy can introduce negative values.
    script:
    if (params.run_eddy) {
        slice_drop_flag=""
        if (params.use_slice_drop_correction) {
            slice_drop_flag="--slice_drop_correction"
        }
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        scil_prepare_eddy_command.py $dwi $bval $bvec $mask\
            --eddy_cmd $params.eddy_cmd --b0_thr $params.b0_thr_extract_b0\
            --encoding_direction $encoding\
            --readout $readout --output_script --fix_seed\
            $slice_drop_flag
        sh eddy.sh
        fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
        mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
    else {
        """
        mv $dwi ${sid}__dwi_corrected.nii.gz
        mv $bvec ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
}

dwi_for_eddy_topup
    .join(gradients_for_eddy_topup)
    .join(topup_files_for_eddy_topup)
    .join(readout_encoding_for_eddy_topup)
    .set{dwi_gradients_mask_topup_files_for_eddy_topup}

process Eddy_Topup {
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0s_corrected),
        file(field), file(movpar), readout, encoding\
        from dwi_gradients_mask_topup_files_for_eddy_topup
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__dwi_corrected.nii.gz", "${sid}__bval_eddy",
        "${sid}__dwi_eddy_corrected.bvec" into\
        dwi_gradients_from_eddy_topup
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_from_eddy_topup
    set sid, "${sid}__bval_eddy", "${sid}__dwi_eddy_corrected.bvec" into\
        gradients_from_eddy_topup
    file "${sid}__b0_bet_mask.nii.gz"

    when:
    rev_b0_count > 0 && params.run_topup

    // Corrected DWI is clipped to ensure there are no negative values
    // introduced by Eddy.
    script:
    if (params.run_eddy) {
        slice_drop_flag=""
        if (params.use_slice_drop_correction)
            slice_drop_flag="--slice_drop_correction"
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        mrconvert $b0s_corrected b0_corrected.nii.gz -coord 3 0 -axes 0,1,2 -nthreads 1
        bet b0_corrected.nii.gz ${sid}__b0_bet.nii.gz -m -R\
            -f $params.bet_topup_before_eddy_f
        scil_prepare_eddy_command.py $dwi $bval $bvec ${sid}__b0_bet_mask.nii.gz\
            --topup $params.prefix_topup --eddy_cmd $params.eddy_cmd\
            --b0_thr $params.b0_thr_extract_b0\
            --encoding_direction $encoding\
            --readout $readout --output_script --fix_seed\
            $slice_drop_flag
        sh eddy.sh
        fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
        mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
    else {
        """
        mv $dwi ${sid}__dwi_corrected.nii.gz
        mv $bvec ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
}

dwi_gradients_from_eddy
    .mix(dwi_gradients_from_eddy_topup)
    .set{dwi_gradients_for_extract_b0}

dwi_from_eddy
    .mix(dwi_from_eddy_topup)
    .set{dwi_for_bet}

gradients_from_eddy
    .mix(gradients_from_eddy_topup)
    .into{gradients_for_resample_b0;
          gradients_for_normalize}

process Extract_B0 {
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradients_for_extract_b0

    output:
    set sid, "${sid}__b0.nii.gz" into b0_for_bet

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    """
}

dwi_for_bet
    .join(b0_for_bet)
    .set{dwi_b0_for_bet}

process Bet_DWI {
    cpus 2

    input:
    set sid, file(dwi), file(b0) from dwi_b0_for_bet

    output:
    set sid, "${sid}__b0_bet.nii.gz", "${sid}__b0_bet_mask.nii.gz" into\
        b0_and_mask_for_crop
    set sid, "${sid}__dwi_bet.nii.gz", "${sid}__b0_bet.nii.gz",
        "${sid}__b0_bet_mask.nii.gz" into dwi_b0_b0_mask_for_n4

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    bet $b0 ${sid}__b0_bet.nii.gz -m -R -f $params.bet_dwi_final_f
    mrcalc $dwi ${sid}__b0_bet_mask.nii.gz -mult ${sid}__dwi_bet.nii.gz -quiet -nthreads 1
    """
}

process N4_DWI {
    cpus 1

    input:
    set sid, file(dwi), file(b0), file(b0_mask)\
        from dwi_b0_b0_mask_for_n4

    output:
    set sid, "${sid}__dwi_n4.nii.gz" into dwi_for_crop

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    N4BiasFieldCorrection -i $b0\
        -o [${sid}__b0_n4.nii.gz, bias_field_b0.nii.gz]\
        -c [300x150x75x50, 1e-6] -v 1
    scil_apply_bias_field_on_dwi.py $dwi bias_field_b0.nii.gz\
        ${sid}__dwi_n4.nii.gz --mask $b0_mask -f
    """
}

dwi_for_crop
    .join(b0_and_mask_for_crop)
    .set{dwi_and_b0_mask_b0_for_crop}

process Crop_DWI {
    cpus 1

    input:
    set sid, file(dwi), file(b0), file(b0_mask) from dwi_and_b0_mask_b0_for_crop

    output:
    set sid, "${sid}__dwi_cropped.nii.gz",
        "${sid}__b0_mask_cropped.nii.gz" into dwi_mask_for_resample
    file "${sid}__b0_cropped.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_crop_volume.py $dwi ${sid}__dwi_cropped.nii.gz -f\
        --output_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0 ${sid}__b0_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0_mask ${sid}__b0_mask_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    """
}

process Resample_DWI {
    cpus 3

    input:
    set sid, file(dwi), file(mask) from dwi_mask_for_resample

    output:
    set sid, "${sid}__dwi_resampled.nii.gz" into\
        dwi_for_resample_b0

    script:
    if (params.run_resample_dwi)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_resample_volume.py $dwi \
            dwi_resample.nii.gz \
            --resolution $params.dwi_resolution \
            --interp  $params.dwi_interpolation
        fslmaths dwi_resample.nii.gz -thr 0 dwi_resample_clipped.nii.gz
        scil_resample_volume.py $mask \
            mask_resample.nii.gz \
            --ref dwi_resample.nii.gz \
            --enforce_dimensions \
            --interp nn
        mrcalc dwi_resample_clipped.nii.gz mask_resample.nii.gz\
            -mult ${sid}__dwi_resampled.nii.gz -quiet -nthreads 1
        """
    else
        """
        mv $dwi ${sid}__dwi_resampled.nii.gz
        """
}

dwi_for_resample_b0
    .join(gradients_for_resample_b0)
    .set{dwi_and_grad_for_resample_b0}

process Resample_B0 {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_and_grad_for_resample_b0

    output:
    file "${sid}__b0_resampled.nii.gz"
    file "${sid}__b0_mask_resampled.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    extract_b0.py $dwi $bval $bvec ${sid}__b0_resampled.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    mrthreshold ${sid}__b0_resampled.nii.gz ${sid}__b0_mask_resampled.nii.gz\
        --abs 0.00001 -nthreads 1
    """
}