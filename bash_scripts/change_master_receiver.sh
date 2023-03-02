#!/bin/bash
REC_ID=$1

echo
echo " changing master receiver to $REC_ID"
echo

rm -f NOISE_TOMOGRAPHY/irec_main_noise
echo "$REC_ID" >> NOISE_TOMOGRAPHY/irec_main_noise 
