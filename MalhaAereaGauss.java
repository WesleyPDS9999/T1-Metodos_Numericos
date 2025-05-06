import java.io.*;
import java.util.*;

public class MalhaAereaGauss {

    public static void main(String[] args) throws IOException {
        String nomeArquivo = "caso200.txt";
        BufferedReader br = new BufferedReader(new FileReader(nomeArquivo));

        // Mapeia o código do aeroporto para o índice da matriz
        Map<String, Integer> indiceAeroporto = new HashMap<>();
        List<String> codigos = new ArrayList<>(); // Lista de códigos de aeroportos
        List<Double> passageirosDiretos = new ArrayList<>(); // Lista de valores de b (passageiros diretos)
        Map<String, Map<String, Double>> fluxo = new HashMap<>(); // fluxo[origem][destino] = porcentagem

        String linha;
        while ((linha = br.readLine()) != null) {
            linha = linha.trim();
            if (linha.isEmpty()) continue;
            String[] partes = linha.split(" ");

            if (partes.length == 2) {
                // Linha com chegada direta (código e número de passageiros)
                String codigo = partes[0];
                double diretos = Double.parseDouble(partes[1]);

                // Adiciona ao mapa e lista se ainda não existir
                if (!indiceAeroporto.containsKey(codigo)) {
                    indiceAeroporto.put(codigo, codigos.size());
                    codigos.add(codigo);
                    passageirosDiretos.add(diretos);
                } else {
                    // Se já existe, atualiza o número de passageiros diretos
                    passageirosDiretos.set(indiceAeroporto.get(codigo), diretos);
                }
            } else if (partes.length == 3) {
                // Linha com fluxo: origem destino porcentagem
                String origem = partes[0];
                String destino = partes[1];
                double perc = Double.parseDouble(partes[2]) / 100.0;

                // Cria entrada no mapa se necessário
                fluxo.putIfAbsent(origem, new HashMap<>());
                fluxo.get(origem).put(destino, perc);

                // Garante que origem e destino existem no índice
                if (!indiceAeroporto.containsKey(origem)) {
                    indiceAeroporto.put(origem, codigos.size());
                    codigos.add(origem);
                    passageirosDiretos.add(0.0); // sem chegada direta, inicialmente
                }
                if (!indiceAeroporto.containsKey(destino)) {
                    indiceAeroporto.put(destino, codigos.size());
                    codigos.add(destino);
                    passageirosDiretos.add(0.0);
                }
            }
        }

        int n = codigos.size(); // número de aeroportos
        double[][] A = new double[n][n]; // matriz de coeficientes
        double[] b = new double[n];      // vetor de constantes (chegadas diretas)

        // Monta o sistema Ax = b
        for (int i = 0; i < n; i++) {
            String destino = codigos.get(i);
            A[i][i] = 1.0;              // começa com x_i (P_i) no lado esquerdo
            b[i] = passageirosDiretos.get(i); // parte direita: passageiros diretos

            // Para cada origem que envia passageiros para o destino
            for (Map.Entry<String, Map<String, Double>> entrada : fluxo.entrySet()) {
                String origem = entrada.getKey();
                Map<String, Double> destinos = entrada.getValue();

                if (destinos.containsKey(destino)) {
                    int j = indiceAeroporto.get(origem);     // posição da origem
                    A[i][j] -= destinos.get(destino);        // subtrai fluxo vindo de origem para destino
                }
            }
        }

        // Resolve o sistema usando eliminação de Gauss
        double[] x = gauss(A, b);

        // Imprime os resultados
        System.out.println("Passageiros em cada aeroporto:");
        double min = Double.MAX_VALUE, max = 0;
        String minCod = "", maxCod = "";

        for (int i = 0; i < n; i++) {
            System.out.printf("%s: %.2f passageiros\n", codigos.get(i), x[i]);

            // Atualiza menor e maior valor
            if (x[i] < min) {
                min = x[i];
                minCod = codigos.get(i);
            }
            if (x[i] > max) {
                max = x[i];
                maxCod = codigos.get(i);
            }
        }

        // Imprime o aeroporto com menor e maior número de passageiros
        System.out.printf("\nMenor população %s: %.2f\n", minCod, min);
        System.out.printf("Maior população %s: %.2f\n", maxCod, max);
    }

    // Função para resolver sistema linear Ax = b com eliminação de Gauss
    public static double[] gauss(double[][] A, double[] b) {
        int n = b.length;

        // Eliminação
        for (int p = 0; p < n; p++) {
            // Escolhe o pivô
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }

            // Troca linhas na matriz e vetor b
            double[] temp = A[p];
            A[p] = A[max];
            A[max] = temp;

            double t = b[p];
            b[p] = b[max];
            b[max] = t;

            // Zera os elementos abaixo do pivô
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // Substituição para trás
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double soma = 0.0;
            for (int j = i + 1; j < n; j++) {
                soma += A[i][j] * x[j];
            }
            x[i] = (b[i] - soma) / A[i][i];
        }

        return x;
    }
}
