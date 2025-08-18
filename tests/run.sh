#!/bin/bash
# run_tests_with_log.sh

LOG_FILE="test_results.log"
ERROR_FILE="test_errors.log"
SUCCESS_COUNT=0
ERROR_COUNT=0
TOTAL_COUNT=0

# Limpar logs anteriores
> "$LOG_FILE"
> "$ERROR_FILE"

echo "Iniciando testes em $(date)" | tee -a "$LOG_FILE"
echo "================================" | tee -a "$LOG_FILE"

while IFS= read -r command; do
    # Pular linhas vazias e comentários
    if [[ -z "$command" || "$command" =~ ^#.* ]]; then
        continue
    fi
    
    TOTAL_COUNT=$((TOTAL_COUNT + 1))
    echo "[$TOTAL_COUNT] Executando: $command" | tee -a "$LOG_FILE"
    
    # Executar comando e capturar saída
    if output=$(eval "$command" 2>&1); then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        echo "✅ SUCESSO" | tee -a "$LOG_FILE"
        echo "$output" >> "$LOG_FILE"
    else
        ERROR_COUNT=$((ERROR_COUNT + 1))
        echo "❌ ERRO" | tee -a "$LOG_FILE" 
        echo "[$TOTAL_COUNT] ERRO: $command" >> "$ERROR_FILE"
        echo "$output" >> "$ERROR_FILE"
        echo "---" >> "$ERROR_FILE"
    fi
    
    echo "------------------------" >> "$LOG_FILE"
done < all_commands.txt

# Resumo final
echo "================================" | tee -a "$LOG_FILE"
echo "Resumo final:" | tee -a "$LOG_FILE"
echo "Total de testes: $TOTAL_COUNT" | tee -a "$LOG_FILE"
echo "Sucessos: $SUCCESS_COUNT" | tee -a "$LOG_FILE"
echo "Erros: $ERROR_COUNT" | tee -a "$LOG_FILE"
echo "Finalizado em $(date)" | tee -a "$LOG_FILE"

if [ $ERROR_COUNT -gt 0 ]; then
    echo "⚠️  Verifique $ERROR_FILE para detalhes dos erros"
fi
