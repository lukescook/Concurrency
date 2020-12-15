//Luke Cook(906771)
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.LinkedList;
import java.util.Queue;

class Group {
    private Queue<Integer> queue = new LinkedList<>();
    private Random random = new Random();
    private int n,m;
    private int i= 0;
    private boolean run = true;

    public void setN(int n) {
        this.n = n;
    }

    public void setM(int m) {
        this.m = m;
    }

    private Integer IncrementReadOdd() {
        if(queue.peek() != null) {
            if (!((queue.peek() % 2) == 0 )) {
                System.out.println("\tOdd number printed");
               return queue.remove();
            }
        }
        return null;
    }

    private Integer IncrementRead(int select) {
        if(queue.peek() != null) {
            if ((queue.peek() % select) == 0 ) {
                System.out.println("\tA number divisible by " + select + " was printed");
                return queue.remove();
            }
        }
        return null;
    }

    private Integer IncrementWrite() {
        int input = random.nextInt();
        queue.add(input);
        System.out.println("\t" + input);
        return input;
    }

    public synchronized Integer Increment(int select) {
        Integer a ;
        if(i>=m && queue.size() ==0){
            run = false;
        }else{
            if(select == 0){
                if(queue.size() < n && i<m){
                    i++;
                    notifyAll();
                    return IncrementWrite();
                    //System.out.print("i " + i + "  m     " + 7);
                } else {
                    try {
                        System.out.println("write wait");
                        notifyAll();
                        wait();
                        System.out.println("write notified");

                    } catch (InterruptedException e) {
                    }
                }
            }else if (queue.size() != 0){
                if(select == 1){
                    a = IncrementReadOdd();
                }else {
                    a = IncrementRead(select);
                }
                if(a != null) {
                    notifyAll();
                    return a;
                } else {
                    try {
                        System.out.println("read wait " + select);
                        wait();
                        System.out.println("read notified " + select);

                    } catch (InterruptedException e) {
                    }
                }
            } else {
                try {
                    System.out.println("read wait " + select);
                    wait();
                    System.out.println("read notified " + select);

                } catch (InterruptedException e) {
                }
            }
        }
        return null;
    }

    public synchronized boolean runState(){
        return run;
    }
}
class ThreadWrite extends Thread {
    private Integer answer;
    private String fileName;
    private int selectValue;
    private Group groupValue;
    private File file;
    private PrintWriter fileOut = null;
    public ThreadWrite(Group setValue, int n, int m) {
        groupValue = setValue;
        groupValue.setN(n);
        groupValue.setM(m);
        this.selectValue = 0;
    }
    public ThreadWrite(Group getValue, int selectValue, String name) {
        fileName = name + ".txt";
        groupValue = getValue;
        this.selectValue = selectValue;
        file = new File(fileName);
        try {
            fileOut = new PrintWriter(file);

        } catch (FileNotFoundException e) {
            System.out.println("Sorry I cannot find " + selectValue);
            System.exit(0);
        }
    }
    public ThreadWrite(Group getValue) {
        groupValue = getValue;
        this.selectValue = 1;
        fileName = "odd-numbers" + ".txt";
        file = new File(fileName);
        try {
            fileOut = new PrintWriter(file);

        } catch (FileNotFoundException e) {
            System.out.println("Sorry I cannot find " + selectValue);
            System.exit(0);
        }
    }

    public void run() {
        while(groupValue.runState()) {
             answer = groupValue.Increment(selectValue);
             if(selectValue != 0 && answer != null){
                 fileOut.println(answer);
             }
        }
        if(selectValue !=0) {
            fileOut.close();
        }
    }
}
public class Count {
    public static void main(String[] Arg) throws InterruptedException {
        Group share = new Group();
        ThreadWrite B = new ThreadWrite(share);
        ThreadWrite C = new ThreadWrite(share, 2,"even-numbers");
        ThreadWrite D = new ThreadWrite(share, 3, "3");
        ThreadWrite A = new ThreadWrite(share,3,10);

        B.start();
        C.start();
        A.start();
        D.start();

        A.join();
        B.join();
        C.join();
        D.join();

        System.out.println("End");


    }
}
