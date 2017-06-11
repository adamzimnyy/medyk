package pwr.adamzimnyy;

import java.io.*;
import java.util.*;

public class Main {

    // liczba wszystkich cech
    static int entrySize = 59;

    final static String[] transposed_files = {
            "D:/MGR/medyk/Zawaly_dane/transpose/inne.txt",
            "D:/MGR/medyk/Zawaly_dane/transpose/ang_prect.txt",
            "D:/MGR/medyk/Zawaly_dane/transpose/ang_prct_2.txt",
            "D:/MGR/medyk/Zawaly_dane/transpose/mi.txt",
            "D:/MGR/medyk/Zawaly_dane/transpose/mi_np.txt"
    };

    final static String transposed_file_all = "D:/MGR/medyk/Zawaly_dane/transpose/all.txt";
    final static String file_all = "D:/MGR/medyk/Zawaly_dane/all.txt";

    final static String[] files = {
            "D:/MGR/medyk/Zawaly_dane/inne.txt",
            "D:/MGR/medyk/Zawaly_dane/ang_prect.txt",
            "D:/MGR/medyk/Zawaly_dane/ang_prct_2.txt",
            "D:/MGR/medyk/Zawaly_dane/mi.txt",
            "D:/MGR/medyk/Zawaly_dane/mi_np.txt"
    };

    final static int METHOD_EUKLIDES = 0;
    final static int METHOD_MANHATTAN = 1;

    final static int NONE = 0;
    final static int STANDARD = 1;
    final static int SCALE = 2;


    static List<DataEntry> dataSet = new ArrayList<>();
    static List<DataEntry> normalDataSet = new ArrayList<>();
    static List<Result> results = new ArrayList<>();
    static int[] priorities = {
            54, 52, 50, 55, 39, 56, 36, 3, 7, 58,
            57, 51, 59, 18, 46, 44, 22, 11, 13, 24,
            47, 53, 4, 26, 43, 23, 20, 48, 31, 8,
            40, 19, 35, 14, 32, 37, 6, 16, 41, 38,
            27, 33, 9, 42, 21, 15, 34, 12, 49, 25,
            1, 5, 29, 17, 45, 2, 28, 10, 30,
    };
    static List<String> traitNames = new ArrayList<>();

    public static void main(String[] args) {


        traitNames.add("age	");
        traitNames.add("sex	");
        traitNames.add("pain_location	");
        traitNames.add("chest_pain_radiation	");
        traitNames.add("pain_character	");
        traitNames.add("onset_of_pain	");
        traitNames.add("number_of_hours_since_onset	");
        traitNames.add("duration_of_the_last_episode	");
        traitNames.add("nausea	");
        traitNames.add("diaphoresis	");
        traitNames.add("palpitations	");
        traitNames.add("dyspnea	");
        traitNames.add("dizziness/syncope	");
        traitNames.add("burping	");
        traitNames.add("palliative_factors	");
        traitNames.add("prior_chest_pain_of_this_type	");
        traitNames.add("physician_consulted_for_prior_pain	");
        traitNames.add("prior_pain_related_to_heart	");
        traitNames.add("prior_pain_due_to_mi	");
        traitNames.add("prior_pain_due_to_angina_prectoris	");
        traitNames.add("prior_mi	");
        traitNames.add("prior_angina_prectoris	");
        traitNames.add("prior_atypical_chest_pain	");
        traitNames.add("congestive_heart_failure	");
        traitNames.add("peripheral_vascular_disease	");
        traitNames.add("hiatal_hernia	");
        traitNames.add("hypertension	");
        traitNames.add("diabetes	");
        traitNames.add("smoker	");
        traitNames.add("diuretics	");
        traitNames.add("nitrates	");
        traitNames.add("beta_blockers	");
        traitNames.add("digitalis	");
        traitNames.add("nonsteroidal_anti-inflammatory	");
        traitNames.add("antacids/h_blockers	");
        traitNames.add("systolic_blood_pressure	");
        traitNames.add("diastolic_blood_pressure	");
        traitNames.add("heart_rate	");
        traitNames.add("respiration_rate	");
        traitNames.add("rales	");
        traitNames.add("cyanosis	");
        traitNames.add("pallor	");
        traitNames.add("systolic_murmur	");
        traitNames.add("diastolic_murmur	");
        traitNames.add("oedema	");
        traitNames.add("_gallop	");
        traitNames.add("s_gallop	");
        traitNames.add("chest_wall_tenderness	");
        traitNames.add("diaphoresis	");
        traitNames.add("new_q_wave	");
        traitNames.add("any_q_wave	");
        traitNames.add("new_st_segment_elevation	");
        traitNames.add("any_st_segment_elevation	");
        traitNames.add("new_st_segment_depression	");
        traitNames.add("any_st_segment_depression	");
        traitNames.add("new_t_wave_inversion	");
        traitNames.add("any_t_wave_inversion	");
        traitNames.add("new_intraventricular_conduction_defect	");
        traitNames.add("any_intraventricular_conduction_defect	");


        for (int i = 0; i < transposed_files.length; i++) {
            File file = new File(transposed_files[i]);
            readFile(file, i + 1);
        }
        System.out.println("File imported. Data size: " + dataSet.size());
        findTraitPriority();

        nearestMeanTest();

        knnTest();

        int[] diagnoses = new int[5];
        for (DataEntry e : dataSet) {
            diagnoses[e.getDiagnosis() - 1]++;
        }
        System.out.println("Diagnoses are: " + Arrays.toString(diagnoses));
        double[] doubles = new double[diagnoses.length];
        for (int i = 0; i < diagnoses.length; i++) {
            doubles[i] = diagnoses[i];
        }
        System.out.println("     standard: " + Arrays.toString(standardizeData(doubles, SCALE)));


        System.out.println("     standard: " + Arrays.toString(standardizeData(new double[]{5, 6, 7, 8, 9, 10}, STANDARD)));

        writeToFile("output.txt");

    }

    static int getNMTableIndex(int normal, int method) {
        if (method == 0 && normal == 0) return 0;
        if (method == 0 && normal == 1) return 1;
        if (method == 0 && normal == 2) return 2;
        if (method == 1 && normal == 0) return 3;
        if (method == 1 && normal == 1) return 4;
        if (method == 1 && normal == 2) return 5;

        throw new ArrayIndexOutOfBoundsException();
    }

    static int getKNNTableIndex(int normal, int method, int k) {
        if (method == 0 && normal == 0 && k == 1) return 0;
        if (method == 0 && normal == 0 && k == 6) return 1;
        if (method == 0 && normal == 0 && k == 15) return 2;
        if (method == 0 && normal == 0 && k == 100) return 3;
        if (method == 1 && normal == 0 && k == 1) return 4;
        if (method == 1 && normal == 0 && k == 6) return 5;
        if (method == 1 && normal == 0 && k == 15) return 6;
        if (method == 1 && normal == 0 && k == 100) return 7;
        if (method == 0 && normal == 1 && k == 1) return 8;
        if (method == 0 && normal == 1 && k == 6) return 9;
        if (method == 0 && normal == 1 && k == 15) return 10;
        if (method == 0 && normal == 1 && k == 100) return 11;
        if (method == 1 && normal == 1 && k == 1) return 12;
        if (method == 1 && normal == 1 && k == 6) return 13;
        if (method == 1 && normal == 1 && k == 15) return 14;
        if (method == 1 && normal == 1 && k == 100) return 15;
        if (method == 0 && normal == 2 && k == 1) return 16;
        if (method == 0 && normal == 2 && k == 6) return 17;
        if (method == 0 && normal == 2 && k == 15) return 18;
        if (method == 0 && normal == 2 && k == 100) return 19;
        if (method == 1 && normal == 2 && k == 1) return 20;
        if (method == 1 && normal == 2 && k == 6) return 21;
        if (method == 1 && normal == 2 && k == 15) return 22;
        if (method == 1 && normal == 2 && k == 100) return 23;


        throw new ArrayIndexOutOfBoundsException("M = " + method + ", N = " + normal + ", K = " + k);
    }

    private static void writeToFile(String filename) {
        try {
            File file = new File(filename);
            if (!file.exists()) file.createNewFile();
            PrintWriter pw = new PrintWriter(new FileWriter(file));
            pw.println("normal = 0, method = METHOD_EUKLIDES");
            double[][] nmTable = new double[entrySize + 1][6];
            double[][] knnTable = new double[entrySize + 1][24];

            for (Result r : results) {
                if (r.k == 0) {
                    nmTable[r.getTraits()][getNMTableIndex(r.getNormalization(), r.getMethod())] = r.getSuccess();
                } else {
                    knnTable[r.getTraits()][getKNNTableIndex(r.getNormalization(), r.getMethod(), r.k)] = r.getSuccess();
                }
            }
            pw.println("\n\n\n ALGORYTM NM \n\n\n");
            StringBuilder sb = new StringBuilder();

            for (int trait = 1; trait <= entrySize; trait++) {
                sb.append(trait).append("\t");
                for (int i = 0; i < 6; i++) {
                    sb.append(nmTable[trait][i]).append("\t");
                }
                sb.append("\n");
            }

            pw.println(sb.toString());
            sb.setLength(0);
            pw.println("\n\n\n ALGORYTM KNN \n\n\n");

            for (int trait = 1; trait <= entrySize; trait++) {
                sb.append(trait).append("\t");
                for (int i = 0; i < 24; i++) {
                    sb.append(knnTable[trait][i]).append("\t");
                }
                sb.append("\n");
            }
            pw.println(sb.toString());


            pw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static void readFile(File file, int diagnosis) {
        int id = 0;
        String readLine = null;
        try {

            BufferedReader b = new BufferedReader(new FileReader(file));
            readLine = "";

            while ((readLine = b.readLine()) != null) {
                if (!readLine.trim().isEmpty()) {
                    DataEntry entry = new DataEntry();
                    entry.setData(readLine.split("\\s"));
                    entry.setDiagnosis(diagnosis);
                    entry.setId(id++);
                    dataSet.add(entry);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        } catch (ArrayIndexOutOfBoundsException ae) {
            System.out.println("ERROR:");
            System.out.println("\t" + Arrays.toString(readLine.split("\\s")));
            System.out.println("\t" + readLine.split("\\s").length);
        }
    }

    public static Map<Integer, DataEntry> calculateMeans(List<DataEntry> data, int maxPriority) {

        Map<Integer, DataEntry> means = new HashMap<>();
        means.clear();
        for (int diagnosis = 1; diagnosis < 6; diagnosis++) {
            DataEntry meanEntry = new DataEntry();
            meanEntry.setDiagnosis(diagnosis);

            double[] avgData = new double[entrySize];
            for (DataEntry e : data) {
                if (e.getDiagnosis() == diagnosis) {
                    avgData = sumArray(avgData, e.getData(), maxPriority);
                }
            }
            avgData = divArray(avgData, entrySize);
            meanEntry.setData(avgData);
            means.put(diagnosis, meanEntry);
        }
        return means;
    }

    public static void nearestMeanTest() {
        System.out.println(" ----- Begin NM test -----");

        int correct = 0;
        int all = 0;
        long startTime = System.currentTimeMillis();
        int predictionStats[] = new int[5];
        int correctPerTrait[] = new int[entrySize];

        for (int normal = 0; normal < 3; normal++) {
            normalDataSet = standardizeData(dataSet, normal);
            System.out.println("normal: " + normal);

            // metoda liczenia odleglosci
            for (int method = 0; method < 2; method++) {
                System.out.println("\tmethod: " + method);

                //liczba cech
                for (int traits = 1; traits <= entrySize; traits++) {
                    if (traits == 1) System.out.print("\t\t");
                    if (traits % 10 == 0) System.out.print(traits);
                    System.out.print(".");
                    if (traits == entrySize) System.out.println();
                    double successRate;
                    correct = 0;
                    all = 0;
                    //5x loop
                    for (int i = 0; i < 5; i++) {

                        // podzial na 2 zestawy
                        List<DataEntry> trainingSample = new ArrayList<>();
                        List<DataEntry> testSample = new ArrayList<>();
                        Collections.shuffle(normal == 1 ? normalDataSet : dataSet);
                        int split = 0;
                        for (DataEntry e : normal == 1 ? normalDataSet : dataSet) {
                            if (split++ % 2 == 0)
                                trainingSample.add(e);
                            else
                                testSample.add(e);
                        }

                        //2cv loop
                        for (int swap = 0; swap < 2; swap++) {

                            // wyliczenie srednich
                            Map<Integer, DataEntry> means = calculateMeans(swap == 0 ? trainingSample : testSample, traits);

                            // przewidywanie diagnozy
                            for (DataEntry current : swap == 0 ? testSample : trainingSample) {
                                double minDistance = Double.MAX_VALUE;
                                int prediction = 0;
                                for (int diagnosis = 1; diagnosis < 6; diagnosis++) {

                                    double distance = getDistance(current.getData(), means.get(diagnosis).getData(), traits, method);
                                    if (distance < minDistance) {
                                        minDistance = distance;
                                        prediction = diagnosis;
                                    }
                                }

                                predictionStats[prediction - 1]++;
                                all++;
                                if (prediction == current.getDiagnosis()) {
                                    correct++;
                                    correctPerTrait[traits - 1]++;
                                }
                            }// prediction loop

                        }//2cv loop

                    } //5x loop
                    //    System.out.println(correct+" ~ "+all);
                    successRate = (double) correct / all;
                    Result r = new Result(traits, method, normal, successRate);
                    results.add(r);

                } // traits loop

            } //method loop

        } // normal loop
        long time = System.currentTimeMillis() - startTime;

        double[] doubles = new double[predictionStats.length];
        for (int i = 0; i < predictionStats.length; i++) {
            doubles[i] = predictionStats[i];
        }

        System.out.println(
                "Test finished." +
                        "\nTesting took " + time + " ms." +
                        " \nTotal results: " + results.size() +
                        " \nPrediction frequency: " + Arrays.toString(predictionStats) +
                        " \n            standard: " + Arrays.toString(standardizeData(doubles, SCALE)) + "\n\n\n"

        );
    }

    static List<Rating> ratings = new ArrayList<>();

    public static void findTraitPriority() {
        int global = 100;
        Map<Integer, Map<Double, Double>> map = new HashMap<>();
        for (int i = 1; i < 6; i++) {
            map.put(i, new HashMap<>());
        }

        // na global mape
        map.put(global, new HashMap<>());

        int[] diagnosisCount = new int[6];
        for (DataEntry e : dataSet) {
            diagnosisCount[e.getDiagnosis()]++;
        }

        for (int trait = 0; trait < entrySize; trait++) {

            for (int i = 1; i < 6; i++) {
                map.get(i).clear();
            }
            map.get(global).clear();
            for (int diagnosis = 1; diagnosis < 6; diagnosis++) {

                // policz ile razy jaka wartość wystąpiła
                for (DataEntry e : dataSet) {
                    if (e.getDiagnosis() == diagnosis) {
                        //dla konkretnej diagnozy
                        if (map.get(diagnosis).containsKey(e.getData()[trait])) {
                            map.get(diagnosis).put(e.getData()[trait], map.get(diagnosis).get(e.getData()[trait]) + 1);
                        } else {
                            map.get(diagnosis).put(e.getData()[trait], 1.0);
                        }

                        //ogólnie
                        if (map.get(global).containsKey(e.getData()[trait])) {
                            map.get(global).put(e.getData()[trait], map.get(global).get(e.getData()[trait]) + 1);
                        } else {
                            map.get(global).put(e.getData()[trait], 1.0);
                        }
                    }
                }

                //policz prawopodobieństwa dla każdej diagnozy
                for (Double key : map.get(diagnosis).keySet()) {
                    map.get(diagnosis).put(key, (map.get(diagnosis).get(key) / diagnosisCount[diagnosis]));
                }
            }

            // policz prawdopodobieństwo całego rozkładu
            for (Double key : map.get(global).keySet()) {
                map.get(global).put(key, map.get(global).get(key) / dataSet.size());
            }


            double rating = 0;

            for (int diagnosis = 1; diagnosis < 6; diagnosis++) {
                double sum = 0;
                //różnice prawdopodobieństwa
                for (Double key : map.get(global).keySet()) {

                    if (map.get(diagnosis).containsKey(key)) {
                        //  System.out.println("\t" + key + ": " + map.get(diagnosis).get(key) + " - " + map.get(global).get(key) + " = " + Math.abs(map.get(diagnosis).get(key) - map.get(global).get(key)));
                        map.get(diagnosis).put(key, Math.abs(map.get(diagnosis).get(key) - map.get(global).get(key)));
                    } else {
                        //   System.out.println("\t" + key + ": " + map.get(diagnosis).get(key) + " - " + map.get(global).get(key));
                        map.get(diagnosis).put(key, Math.abs(map.get(global).get(key)));
                    }
                    sum += map.get(diagnosis).get(key);
                }
                rating += sum * (double) diagnosisCount[diagnosis] / dataSet.size();

            }
            Rating r = new Rating();
            r.rating = rating;
            r.trait = trait;
            r.name = traitNames.get(trait);
            ratings.add(r);
        }

        // posortuj po ocenie, przypisz priorytety po kolei
        ratings.sort((o1, o2) -> new Double(o2.rating).compareTo(o1.rating));
        for (int prio = 0; prio < entrySize; prio++) {
            ratings.get(prio).priority = prio;
        }

        // posortuj w kolejnosc domysla, spisz priorytety
        ratings.sort((o1, o2) -> new Integer(o1.trait).compareTo(o2.trait));

        for (int i = 0; i < entrySize; i++) {
            priorities[i] = ratings.get(i).priority;
        }
    }

    public static double[] sumArray(double[] a, double[] b, int traitPriority) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            //  if (priorities[i] <= traitPriority)
            c[i] = a[i] + b[i];
        }
        return c;
    }

    public static double[] subArray(double[] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            c[i] = a[i] - b[i];
        }
        return c;
    }

    public static double[] divArray(double[] a, int div) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            c[i] = a[i] / div;
        }
        return c;
    }

    public static double[] divArray(int[] a, int div) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            c[i] = (double) a[i] / div;
        }
        return c;
    }

    public static double getDistance(double[] a, double[] b, int traitPriority, int method) {
        double distance = 0;
        if (method == METHOD_EUKLIDES) {
            double diff[] = subArray(a, b);
            for (int i = 0; i < diff.length; i++) {
                if (priorities[i] <= traitPriority)
                    distance += Math.pow(diff[i], 2);
            }
            distance = Math.sqrt(distance);

        } else if (method == METHOD_MANHATTAN) {
            double sum[] = subArray(a, b);
            for (int i = 0; i < sum.length; i++) {
                if (priorities[i] <= traitPriority)
                    distance += Math.abs(sum[i]);
            }
        }
        return distance;
    }


    public static double getDistance(double[] a, double[] b, int method) {
        double distance = 0;
        if (method == METHOD_EUKLIDES) {
            double diff[] = subArray(a, b);
            for (int i = 0; i < diff.length; i++) {
                distance += Math.pow(diff[i], 2);
            }
            distance = Math.sqrt(distance);

        } else if (method == METHOD_MANHATTAN) {
            double sum[] = subArray(a, b);
            for (int i = 0; i < sum.length; i++) {
                distance += Math.abs(sum[i]);
            }
        }
        return distance;
    }


    public static double[] standardizeData(double[] data, int mode) {
        double[] result = new double[data.length];
        System.arraycopy(data, 0, result, 0, data.length);

        if (mode == STANDARD) {
            double mean, std, sum = 0, temp = 0;
            for (double a : data)
                sum += a;
            mean = sum / data.length;

            for (double a : data)
                temp += (a - mean) * (a - mean);
            std = Math.sqrt(temp / (data.length - 1));

            for (int i = 0; i < result.length; i++) {
                result[i] = (result[i] - mean) / std;
            }

        } else if (mode == SCALE) {
            double min = Double.MAX_VALUE;
            double max = 0;
            for (int j = 0; j < data.length; j++) {
                max = max < result[j] ? result[j] : max;
                min = min > result[j] ? result[j] : min;
            }
            for (int j = 0; j < data.length; j++) {
                result[j] = (result[j] - min) / (max - min);
            }
        }
        return result;
    }


    public static List<DataEntry> standardizeData(List<DataEntry> data, int mode) {
        System.out.println("Standardize data in mode " + mode);

        List<DataEntry> result = new ArrayList<>(data);

        if (mode == STANDARD) {
            double[] mean = new double[entrySize];
            double[] std = new double[entrySize];

            for (int i = 0; i < data.size(); i++) {
                for (int j = 0; j < entrySize; j++) {
                    mean[j] += data.get(i).getData()[j];
                    if (i == data.size() - 1) {
                        mean[j] /= data.size();

                    }
                }
            }

            for (int j = 0; j < entrySize; j++) {
                for (int i = 0; i < data.size(); i++) {
                    std[j] += Math.pow(data.get(i).getData()[j] - mean[j], 2);
                    if (i == data.size() - 1) {
                        std[j] = Math.sqrt(std[j] / data.size());
                    }
                }
            }

            for (DataEntry e : result) {
                for (int i = 0; i < e.getData().length; i++) {
                    e.getData()[i] = (e.getData()[i] - mean[i]) / std[i];
                }
            }
        } else if (mode == SCALE) {
            double min[] = new double[entrySize];
            double max[] = new double[entrySize];
            for (int j = 0; j < entrySize; j++) {
                min[j] = 0;
                max[j] = 0;
            }
            for (int j = 0; j < entrySize; j++) {
                for (DataEntry e : result) {
                    if (max[j] < e.getData()[j]) {
                        max[j] = e.getData()[j];
                    }
                   /* if (min[j] > e.getData()[j]) {
                        min[j] = e.getData()[j];
                    }*/
                }
            }
            for (int j = 0; j < entrySize; j++) {
                for (DataEntry e : result) {
                    e.getData()[j] = (e.getData()[j] - min[j]) / (max[j] - min[j]);
                }
            }
        }
        return result;
    }


    //region ALGORYTM K NAJBLIZSZYCH SASIADOW

    public static void knnTest() {
        System.out.println(" ----- Begin KNN test -----");

        int[] kValues = {1, 6, 15, 100};
        int correct = 0;
        int all = 0;
        long startTime = System.currentTimeMillis();
        int predictionStats[] = new int[5];

        for (int normal = 0; normal < 3; normal++) {
            normalDataSet = standardizeData(dataSet, normal);
            System.out.println("normal: " + normal);
            // metoda liczenia odleglosci
            for (int method = 0; method < 2; method++) {
                System.out.println("\tmethod: " + method);

                //liczba cech
                for (int traits = 1; traits <= entrySize; traits++) {
                    if (traits == 1) System.out.print("\t\t");
                    if (traits % 10 == 0) System.out.print(traits);
                    System.out.print(".");
                    if (traits == entrySize) System.out.println();
                    double successRate;
                    for (int k = 0; k < 4; k++) {
                        correct = 0;
                        all = 0;
                        //5x loop
                        for (int i = 0; i < 5; i++) {

                            // podzial na 2 zestawy
                            List<DataEntry> trainingSample = new ArrayList<>();
                            List<DataEntry> testSample = new ArrayList<>();
                            Collections.shuffle(normal == 1 ? normalDataSet : dataSet);
                            int split = 0;
                            for (DataEntry e : normal == 1 ? normalDataSet : dataSet) {
                                if (split++ % 2 == 0)
                                    trainingSample.add(e);
                                else
                                    testSample.add(e);
                            }

                            //2cv loop
                            for (int swap = 0; swap < 2; swap++) {

                                for (DataEntry current : swap == 1 ? trainingSample : testSample) {

                                    int[] diagnoses = new int[6];
                                    List<DataEntry> neighbours = getNeighbours(current,
                                            swap == 1 ? testSample : trainingSample, method, traits, kValues[k]);
                                    for (DataEntry nei : neighbours) {
                                        diagnoses[nei.getDiagnosis()]++;
                                    }
                                    int max = 0;
                                    int prediction = 0;
                                    for (int diag = 1; diag < 6; diag++) {
                                        if (diagnoses[diag] > max) {
                                            max = diagnoses[diag];
                                            prediction = diag;
                                        }
                                    }
                                    if (prediction == current.getDiagnosis()) {
                                        correct++;
                                    }
                                    all++;
                                    predictionStats[prediction - 1]++;
                                }
                            }//2cv loop

                        } //5x loop
                        //    System.out.println(correct+" ~ "+all);
                        successRate = (double) correct / all;
                        Result r = new Result(traits, method, normal, successRate, kValues[k]);
                        results.add(r);

                    } // k loop
                } // traits loop

            } //method loop

        } // normal loop
        long time = System.currentTimeMillis() - startTime;


        double[] doubles = new double[predictionStats.length];
        for (int i = 0; i < predictionStats.length; i++) {
            doubles[i] = predictionStats[i];
        }


        System.out.println(
                "Test finished." +
                        "\nTesting took " + time + " ms." +
                        " \nTotal results: " + results.size() +
                        " \nPrediction frequency: " + Arrays.toString(predictionStats) +
                        " \n            standard: " + Arrays.toString(standardizeData(doubles, SCALE)) + "\n\n\n"

        );
    }

    public static List<DataEntry> getNeighbours(DataEntry testSample, List<DataEntry> trainingSet, int method, int traits, int k) {
        List<DataEntry> neighbours = new ArrayList<>();

        List<Distance> distances = new ArrayList<>();
        for (DataEntry e : trainingSet) {
            double distance = getDistance(e.getData(), testSample.getData(), traits, method);
            distances.add(new Distance(distance, e));
        }
        distances.sort((o1, o2) -> new Double(o1.distance).compareTo(o2.distance));
        for (int i = 0; i < k; i++) {
            neighbours.add(distances.get(i).getPoint());
        }
        return neighbours;
    }
}




