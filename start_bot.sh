#!/bin/bash
echo "Starting Gap Fade Trading Bot..."
echo "Opening TWS first - make sure you log in"
echo "Bot will trade at 14:35 UK time"
echo ""
cd ~/Desktop/Evevrything/Projects/gapfade
caffeinate -i python3 gap_fade_trader.py
